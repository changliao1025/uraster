import os
import logging
import traceback
from typing import Optional, Tuple, List, Dict, Any, Union
import numpy as np
from osgeo import gdal, ogr
from uraster.classes.sraster import sraster
gdal.UseExceptions()
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.international_date_line_utility import split_international_date_line_polygon_coordinates, check_cross_international_date_line_polygon
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
from pyearth.gis.geometry.extract_unique_vertices_and_connectivity import extract_unique_vertices_and_connectivity
# Try to import psutil for memory monitoring (optional)
try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False

# Set up logging
logger = logging.getLogger(__name__)
crs = "EPSG:4326"
# Utility functions for common operations

def _log_memory_usage(stage: str, iFlag_verbose_in: bool = False) -> None:
    """
    Log current memory usage if psutil is available.

    Args:
        stage: Description of the current processing stage
        iFlag_verbose_in: Whether to log memory information
    """
    if not PSUTIL_AVAILABLE or not iFlag_verbose_in:
        return

    try:
        process = psutil.Process()
        memory_info = process.memory_info()
        memory_mb = memory_info.rss / (1024 * 1024)
        logger.info(f"Memory usage at {stage}: {memory_mb:.1f} MB")
    except Exception as e:
        logger.debug(f"Could not get memory usage: {e}")

def check_geometry_validity(sFilename_source_mesh: str, iFlag_verbose_in: bool = False) -> bool:
    """
    Comprehensive check of all geometries in a mesh vector file.

    Consolidates all polygon geometry validation including:
    - Coordinate range validation (-180 to 180 for lon, -90 to 90 for lat)
    - OGR geometry validity checks
    - International Date Line crossing detection
    - Multipolygon part validation
    - Minimum vertex count checks

    Args:
        sFilename_source_mesh (str): Path to the source mesh vector file
        iFlag_verbose_in (bool): If True, print detailed progress messages

    Returns:
        bool: True if all geometries are valid, False otherwise
    """
    try:
        pDataset = ogr.Open(sFilename_source_mesh, 0)  # Read-only
        if pDataset is None:
            logger.error(f'Failed to open file: {sFilename_source_mesh}')
            return False

        if iFlag_verbose_in:
            logger.info(f'Successfully opened mesh file: {sFilename_source_mesh}')

        # Get the first layer
        pLayer = pDataset.GetLayer(0)
        if pLayer is None:
            logger.error('Failed to get layer from the dataset.')
            pDataset = None
            return False

        # Get layer information
        pLayerDefn = pLayer.GetLayerDefn()
        if pLayerDefn is None:
            logger.error('Failed to get layer definition.')
            pDataset = None
            return False

        nFeatures = pLayer.GetFeatureCount()
        if nFeatures == 0:
            logger.warning('Layer contains no features.')
            pDataset = None
            return False

        if iFlag_verbose_in:
            logger.info(f'Validating geometries for {nFeatures} features...')

        # Process features with comprehensive validation
        pLayer.ResetReading()
        iFeature_index = 0
        invalid_geometry_count = 0
        valid_geometry_count = 0

        for pFeature in pLayer:
            if pFeature is None:
                invalid_geometry_count += 1
                iFeature_index += 1
                continue

            pGeometry = pFeature.GetGeometryRef()
            if pGeometry is None:
                logger.warning(f'Feature {iFeature_index} has no geometry')
                invalid_geometry_count += 1
                iFeature_index += 1
                continue

            # Skip GDAL geometry validation as it cannot handle IDL-crossing cells

            sGeometry_type = pGeometry.GetGeometryName()

            aCoord = get_geometry_coordinates(pGeometry)
            #also check longitude and latitude range -180-180
            if aCoord is None or len(aCoord) < 3:
                logger.warning(f'Feature {iFeature_index}: Invalid or insufficient coordinates')
                invalid_geometry_count += 1
                iFeature_index += 1
                continue
            if (np.any(aCoord[:,0] < -180) or np.any(aCoord[:,0] > 180) or
                np.any(aCoord[:,1] < -90) or np.any(aCoord[:,1] > 90)):
                logger.warning(f'Feature {iFeature_index}: Coordinates out of valid range')
                invalid_geometry_count += 1
                iFeature_index += 1
                continue

            if sGeometry_type == 'POLYGON':
                if not _validate_polygon_geometry(pGeometry, iFeature_index, iFlag_verbose_in):
                    invalid_geometry_count += 1
                else:
                    valid_geometry_count += 1

            elif sGeometry_type == 'MULTIPOLYGON':
                if not _validate_multipolygon_geometry(pGeometry, iFeature_index, iFlag_verbose_in):
                    invalid_geometry_count += 1
                else:
                    valid_geometry_count += 1

            elif sGeometry_type in ['POINT', 'LINESTRING']:
                logger.warning(f'Feature {iFeature_index}: Geometry type {sGeometry_type} not supported for mesh processing')
                invalid_geometry_count += 1

            else:
                logger.warning(f'Feature {iFeature_index}: Unknown geometry type {sGeometry_type}')
                invalid_geometry_count += 1

            iFeature_index += 1

        # Cleanup
        pDataset = None

        # Report validation results
        total_features = valid_geometry_count + invalid_geometry_count
        success_rate = (valid_geometry_count / total_features * 100) if total_features > 0 else 0

        if iFlag_verbose_in or invalid_geometry_count > 0:
            logger.info(f'Geometry validation summary:')
            logger.info(f'  - Total features processed: {total_features}')
            logger.info(f'  - Valid geometries: {valid_geometry_count}')
            logger.info(f'  - Invalid geometries: {invalid_geometry_count}')
            logger.info(f'  - Success rate: {success_rate:.1f}%')

        if invalid_geometry_count > 0:
            logger.warning('Found invalid geometries. The program will attempt to fix them.')
            return False

        if iFlag_verbose_in:
            logger.info('All geometries passed validation')

        return True

    except Exception as e:
        logger.error(f'Error in check_geometry_validity: {str(e)}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return False

def check_mesh_quality(sFilename_mesh_in, iFlag_verbose_in=False):
    if not check_geometry_validity(sFilename_mesh_in, iFlag_verbose_in=iFlag_verbose_in):
        # we need to fix the mesh using the IDL splitting utility
        # Make the filename adjustment more flexible to handle any format
        # Get the file extension and base name
        file_base, file_ext = os.path.splitext(sFilename_mesh_in)
        sFilename_source_mesh_fixed = f"{file_base}_fixed{file_ext}"
        fix_mesh_longitude_range_and_idl_crossing(sFilename_mesh_in, sFilename_source_mesh_fixed)
        return sFilename_source_mesh_fixed
    return sFilename_mesh_in

def _validate_polygon_geometry(pGeometry, feature_id, iFlag_verbose_in=False):
    """
    Validate a single polygon geometry including coordinate range and IDL checks.

    Args:
        pGeometry: OGR Polygon geometry
        feature_id: Feature identifier for logging
        iFlag_verbose_in: Verbose logging flag

    Returns:
        bool: True if geometry is valid, False otherwise
    """
    try:
        # Get coordinates
        aCoord = get_geometry_coordinates(pGeometry)
        if aCoord is None or len(aCoord) < 3:
            logger.warning(f'Feature {feature_id}: Invalid or insufficient coordinates for polygon')
            return False

        # Validate coordinate bounds
        lons = aCoord[:, 0]
        lats = aCoord[:, 1]

        # Check coordinate ranges
        if (np.any(lons < -180) or np.any(lons > 180) or
            np.any(lats < -90) or np.any(lats > 90)):
            logger.warning(f'Feature {feature_id}: Coordinates outside valid range')
            logger.warning(f'  Longitude range: {lons.min():.3f} to {lons.max():.3f}')
            logger.warning(f'  Latitude range: {lats.min():.3f} to {lats.max():.3f}')
            return False

        # Check for International Date Line crossing (this is allowed but logged)
        iCross_idl, dummy = check_cross_international_date_line_polygon(aCoord)
        if iCross_idl:
            if iFlag_verbose_in:
                logger.info(f'Feature {feature_id}: Polygon crosses International Date Line (valid)')
        else:
            #use gdal geometry validity check only when it does not cross IDL
            if not pGeometry.IsValid():
                logger.warning(f'Feature {feature_id}: Polygon geometry is invalid according to OGR')
                return False
        return True

    except Exception as e:
        logger.warning(f'Feature {feature_id}: Error validating polygon geometry: {str(e)}')
        return False

def _validate_multipolygon_geometry(pGeometry, feature_id, iFlag_verbose_in=False):
    """
    Validate a multipolygon geometry by checking all constituent polygons.

    Args:
        pGeometry: OGR MultiPolygon geometry
        feature_id: Feature identifier for logging
        iFlag_verbose_in: Verbose logging flag

    Returns:
        bool: True if all parts are valid, False otherwise
    """
    try:
        valid_parts = 0
        total_parts = pGeometry.GetGeometryCount()

        if total_parts == 0:
            logger.warning(f'Feature {feature_id}: Multipolygon has no parts')
            return False

        for iPart in range(total_parts):
            pPolygon_part = pGeometry.GetGeometryRef(iPart)
            if pPolygon_part is None:
                logger.warning(f'Feature {feature_id}: Multipolygon part {iPart} is None')
                continue

            # Skip GDAL geometry validation as it cannot handle IDL-crossing cells

            # Validate coordinates of this part
            aCoord_part = get_geometry_coordinates(pPolygon_part)
            if aCoord_part is None or len(aCoord_part) < 3:
                logger.warning(f'Feature {feature_id}: Multipolygon part {iPart} has insufficient coordinates')
                continue

            # Check coordinate bounds for this part
            lons_part = aCoord_part[:, 0]
            lats_part = aCoord_part[:, 1]

            if (np.any(lons_part < -180) or np.any(lons_part > 180) or
                np.any(lats_part < -90) or np.any(lats_part > 90)):
                logger.warning(f'Feature {feature_id}: Multipolygon part {iPart} has coordinates outside valid range')
                continue

            valid_parts += 1

        if valid_parts == 0:
            logger.warning(f'Feature {feature_id}: No valid parts found in multipolygon')
            return False

        if iFlag_verbose_in and valid_parts < total_parts:
            logger.info(f'Feature {feature_id}: Multipolygon has {valid_parts}/{total_parts} valid parts')

        return True

    except Exception as e:
        logger.warning(f'Feature {feature_id}: Error validating multipolygon geometry: {str(e)}')
        return False

def get_polygon_list(
    sFilename_source_mesh: str,
    iFlag_verbose_in: bool = False,
    sField_unique_id: str = 'cellid'
) -> Optional[Tuple[List[Tuple[Union[int, str], str]], List[float], Optional[str]]]:
    """
    Extract polygon geometries and areas from mesh vector file.

    Processes mesh features, handles International Date Line (IDL) crossing polygons,
    and returns polygon WKT strings with associated areas and projection information.

    Args:
        sFilename_source_mesh (str): Path to the source mesh vector file
        iFlag_verbose_in (bool, optional): If True, print detailed progress messages.
            Default is False.
        sField_unique_id (str, optional): Name of the field containing unique cell IDs.
            Default is 'cellid'.

    Returns:
        Optional[Tuple[List[Tuple[Union[int, str], str]], List[float], Optional[str]]]:
            - List of (cellid, wkt_string) tuples for each polygon
            - List of polygon areas in square degrees
            - Source projection WKT string
            Returns None on failure.

    Raises:
        ValueError: If input parameters are invalid
    """
    # Input validation
    if not isinstance(sFilename_source_mesh, str) or not sFilename_source_mesh.strip():
        logger.error("Invalid mesh filename provided")
        return None

    if not os.path.exists(sFilename_source_mesh):
        logger.error(f"Mesh file does not exist: {sFilename_source_mesh}")
        return None

    if iFlag_verbose_in:
        logger.info(
            "get_polygon_list: Pre-fetching features and analyzing geometries...")

    aPolygon = []
    aArea = []
    pDataset_mesh = None
    pLayer_mesh = None
    pSpatialRef_source = None

    try:
        # Open the mesh vector file
        pDataset_mesh = ogr.Open(sFilename_source_mesh, 0)  # 0 means read-only
        if pDataset_mesh is None:
            logger.error(
                f"Failed to open mesh dataset: {sFilename_source_mesh}")
            return None

        pLayer_mesh = pDataset_mesh.GetLayer(0)
        if pLayer_mesh is None:
            logger.error("Failed to get layer from mesh dataset")
            return None

        nFeature = pLayer_mesh.GetFeatureCount()
        if nFeature <= 0:
            logger.warning("No features found in mesh dataset")
            return [], [], None

        if iFlag_verbose_in:
            logger.info(f"Found {nFeature} features in mesh dataset")

        pSpatialRef_source = pLayer_mesh.GetSpatialRef()
        sProjection_source_wkt = pSpatialRef_source.ExportToWkt() if pSpatialRef_source else None

        if sProjection_source_wkt is None:
            logger.warning("No spatial reference found in mesh dataset")

        # Process features
        pLayer_mesh.ResetReading()
        i = 0
        processed_count = 0
        error_count = 0

        for pFeature_mesh in pLayer_mesh:
            if pFeature_mesh is None:
                error_count += 1
                continue

            try:
                # Get geometry (validation already done by check_geometry_validity)
                pPolygon = pFeature_mesh.GetGeometryRef()
                if pPolygon is None:
                    error_count += 1
                    i += 1
                    continue

                sGeometry_type = pPolygon.GetGeometryName()
                # Read cellid from current feature with error handling
                try:
                    # Handle both string and integer field types
                    pField_defn = pLayer_mesh.GetLayerDefn().GetFieldDefn(
                        pLayer_mesh.GetLayerDefn().GetFieldIndex(sField_unique_id))
                    if pField_defn.GetType() == ogr.OFTString:
                        current_cellid = pFeature_mesh.GetFieldAsString(sField_unique_id)
                    else:
                        current_cellid = pFeature_mesh.GetFieldAsInteger(sField_unique_id)

                    if current_cellid is None or current_cellid == '':
                        current_cellid = i  # Use feature index as fallback
                except Exception as field_error:
                    logger.warning(
                        f"Error reading {sField_unique_id} for feature {i}: {field_error}")
                    current_cellid = i

                if sGeometry_type == "POLYGON":
                    try:
                        aCoord = get_geometry_coordinates(pPolygon)
                        if aCoord is None or len(aCoord) < 3:
                            error_count += 1
                            i += 1
                            continue

                        #no more need to split IDL crossing polygon here, as we have handled it in check_geometry_validity
                        # Regular polygon (no IDL crossing)
                        try:
                            dArea = calculate_polygon_area(
                                aCoord[:, 0], aCoord[:, 1])
                            wkt = pPolygon.ExportToWkt()
                            aPolygon.append((current_cellid, wkt))
                            aArea.append(dArea)
                            processed_count += 1
                        except Exception as area_error:
                            logger.warning(
                                f"Error calculating area for feature {i}: {area_error}")
                            error_count += 1

                    except Exception as polygon_error:
                        logger.warning(
                            f"Error processing polygon feature {i}: {polygon_error}")
                        error_count += 1

                elif sGeometry_type == "MULTIPOLYGON":
                    try:
                        dArea = 0.0
                        for iPart in range(pPolygon.GetGeometryCount()):
                            pPolygon_part = pPolygon.GetGeometryRef(iPart)
                            if pPolygon_part is None:
                                continue
                            try:
                                aCoords_part = get_geometry_coordinates(
                                    pPolygon_part)
                                if aCoords_part is not None and len(aCoords_part) >= 3:
                                    dArea += calculate_polygon_area(
                                        aCoords_part[:, 0], aCoords_part[:, 1])
                            except Exception as part_error:
                                logger.warning(
                                    f"Error processing multipolygon part {iPart} of feature {i}: {part_error}")
                                continue

                        wkt = pPolygon.ExportToWkt()
                        aPolygon.append((current_cellid, wkt))
                        aArea.append(dArea)
                        processed_count += 1

                    except Exception as multipolygon_error:
                        logger.warning(
                            f"Error processing multipolygon feature {i}: {multipolygon_error}")
                        error_count += 1
                else:
                    logger.warning(
                        f"Unsupported geometry type '{sGeometry_type}' for feature {i}")
                    error_count += 1

            except Exception as feature_error:
                logger.warning(
                    f"Error processing feature {i}: {feature_error}")
                error_count += 1

            i += 1

            # Progress reporting during feature pre-processing
            if i % 1000 == 0 and iFlag_verbose_in:
                logger.info(
                    f"Pre-processed {i} features... ({processed_count} successful, {error_count} errors)")

        # Final summary
        if iFlag_verbose_in:
            logger.info(f"get_polygon_list: Pre-processing completed.")
            logger.info(f"  Total features processed: {i}")
            logger.info(f"  Successfully processed: {processed_count}")
            logger.info(f"  Errors/skipped: {error_count}")
            logger.info(
                f"  Success rate: {(processed_count/i*100):.1f}%" if i > 0 else "  Success rate: 0%")

        return aPolygon, aArea, sProjection_source_wkt

    except Exception as e:
        logger.error(f"Error in get_polygon_list: {str(e)}")
        logger.error(f"Traceback: {traceback.format_exc()}")
        return None

    finally:
        # Cleanup resources
        try:
            if pSpatialRef_source is not None:
                pSpatialRef_source = None
        except Exception as e:
            logger.warning(f"Error cleaning up spatial reference: {e}")

        try:
            if pLayer_mesh is not None:
                pLayer_mesh = None
        except Exception as e:
            logger.warning(f"Error cleaning up layer: {e}")

        try:
            if pDataset_mesh is not None:
                pDataset_mesh = None
        except Exception as e:
            logger.warning(f"Error cleaning up dataset: {e}")

def get_unique_values_from_rasters(aFilename_raster: str,
                                    dMissing_value: float,
                                    band_index: int = 1,
                                    iFlag_verbose_in: bool = False) -> Optional[List[float]]:
    """
    Extract unique values from a raster band.

    Args:
        sFilename_raster (str): Path to the raster file
        band_index (int, optional): Band index to read (1-based). Default is

    """
    aUnique_values = set()
    for sFilename in aFilename_raster:
        pRaster = sraster(sFilename)
        if pRaster is not None:
            unique_values = pRaster.get_unique_values(band_index, dMissing_value, iFlag_verbose_in)
            if unique_values is not None:
                aUnique_values.update(unique_values)

    return list(aUnique_values) if aUnique_values else None

def rebuild_mesh_topology(sFilename_mesh_in, iFlag_verbose_in=False, sField_unique_id=None):
    """
    Rebuild mesh topology from source mesh file by extracting vertices,
    connectivity, and centroids for unstructured mesh processing.

    Args:
        sFilename_mesh_in (str): Path to the source mesh file (GeoJSON, Shapefile, etc.)
        iFlag_verbose_in (bool, optional): If True, print detailed progress messages.
            If False, only print error messages. Default is False.
        sField_unique_id (str, optional): Field name for unique cell IDs. If None, uses first field or feature index.
            Note: Field is always treated as integer type since setup_mesh_cellid() enforces this.

    Returns:
        dict: Comprehensive mesh topology information with keys:
            - 'vertices_longitude': np.ndarray of unique vertex longitudes
            - 'vertices_latitude': np.ndarray of unique vertex latitudes
            - 'connectivity': np.ndarray connectivity matrix
            - 'cell_centroids_longitude': np.ndarray of cell centroid longitudes
            - 'cell_centroids_latitude': np.ndarray of cell centroid latitudes
            - 'cell_ids': np.ndarray of cell IDs
            - 'area_min': float minimum cell area
            - 'area_max': float maximum cell area
            - 'area_mean': float mean cell area
            - 'max_vertices_per_cell': int maximum vertices per cell
            - 'num_cells': int total number of cells
            - 'num_vertices': int total number of unique vertices
            - 'success': bool whether processing was successful
        Returns None on failure.
    """
    try:
        # Open the input data source
        pDataset = ogr.Open(sFilename_mesh_in, 0)  # Read-only
        if pDataset is None:
            logger.error(f'Failed to open file: {sFilename_mesh_in}')
            return None
        if iFlag_verbose_in:
            logger.info(f'Successfully opened mesh file: {sFilename_mesh_in}')
        # Get the first layer
        pLayer = pDataset.GetLayer(0)
        if pLayer is None:
            logger.error('Failed to get layer from the dataset.')
            pDataset = None
            return None
        # Get layer information
        pLayerDefn = pLayer.GetLayerDefn()
        if pLayerDefn is None:
            logger.error('Failed to get layer definition.')
            pDataset = None
            return None
        nFeatures = pLayer.GetFeatureCount()

        iFieldCount = pLayerDefn.GetFieldCount()
        if nFeatures == 0:
            logger.warning('Layer contains no features.')
            pDataset = None
            return None
        aCellID= []  # Will be populated dynamically as features are processed
        if iFlag_verbose_in:
            logger.info(f'Processing {nFeatures} features with {iFieldCount} fields')
        # Get the first field name (assuming it contains the data variable)
        if sField_unique_id is None:
            sVariable = pLayerDefn.GetFieldDefn(0).GetName() if iFieldCount > 0 else None
        else:
            sVariable = sField_unique_id
        # Initialize lists for storing geometry data
        lons_list = []
        lats_list = []
        aArea_list = []
        # Process features with enhanced error handling
        pLayer.ResetReading()
        iFeature_index = 0
        invalid_geometry_count = 0
        iCount_multipolygon_cells = 0
        for pFeature in pLayer:
            if pFeature is None:
                continue
            pGeometry = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry.GetGeometryName()
            if sGeometry_type == 'POLYGON':
                try:
                    # Get coordinates of the polygon (validation already done by check_geometry_validity)
                    aCoord = get_geometry_coordinates(pGeometry)
                    if aCoord is not None and len(aCoord) >= 3:
                        lons = aCoord[:, 0]
                        lats = aCoord[:, 1]
                        lons_list.append(lons)
                        lats_list.append(lats)
                        # Calculate polygon area
                        try:
                            dArea = calculate_polygon_area(lons, lats)
                            aArea_list.append(dArea)
                        except Exception as area_error:
                            logger.warning(f'Could not calculate area: {area_error}')
                            aArea_list.append(0.0)
                        # Get field data (always integer since setup_mesh_cellid enforces it)
                        if sVariable:
                            try:
                                field_value = pFeature.GetFieldAsInteger(sVariable)
                                aCellID.append(int(field_value))
                            except (ValueError, TypeError, AttributeError) as e:
                                logger.warning(f'Could not read integer field value: {e}')
                                aCellID.append(len(aCellID))
                        else:
                            aCellID.append(len(aCellID))
                    else:
                        invalid_geometry_count += 1
                except Exception as e:
                    logger.warning(f'Error processing feature: {str(e)}')
                    invalid_geometry_count += 1
            elif sGeometry_type == 'MULTIPOLYGON':
                try:
                    # Process multipolygon by extracting each constituent polygon
                    if iFlag_verbose_in:
                        logger.info('Processing multipolygon feature')
                    multipolygon_processed = False
                    iCount_multipolygon_cells += 1
                    for iPart in range(pGeometry.GetGeometryCount()):
                        pPolygon_part = pGeometry.GetGeometryRef(iPart)
                        if pPolygon_part is None:
                            logger.warning('Multipolygon part is None')
                            continue
                        # Get coordinates of the polygon part (validation already done by check_geometry_validity)
                        aCoord_part = get_geometry_coordinates(pPolygon_part)
                        if aCoord_part is not None and len(aCoord_part) >= 3:
                            lons_part = aCoord_part[:, 0]
                            lats_part = aCoord_part[:, 1]
                            lons_list.append(lons_part)
                            lats_list.append(lats_part)
                            # Calculate polygon area for this part
                            try:
                                dArea_part = calculate_polygon_area(lons_part, lats_part)
                                aArea_list.append(dArea_part)
                            except Exception as area_error:
                                logger.warning(f'Could not calculate area for multipolygon part: {area_error}')
                                aArea_list.append(0.0)
                            multipolygon_processed = True
                        else:
                            logger.warning('Failed to extract coordinates from multipolygon part')
                    if not multipolygon_processed:
                        logger.warning('No valid parts found in multipolygon feature')
                        invalid_geometry_count += 1
                except Exception as e:
                    logger.warning(f'Error processing multipolygon feature: {str(e)}')
                    invalid_geometry_count += 1
            elif sGeometry_type in ['POINT', 'LINESTRING']:
                logger.warning(f'Geometry type {sGeometry_type} not supported in feature {iFeature_index}, skipping')
                invalid_geometry_count += 1
            else:
                logger.warning(f'Unknown geometry type {sGeometry_type} in feature {iFeature_index}, skipping')
                invalid_geometry_count += 1
        # Report processing statistics
        valid_mesh_cells = len(lons_list)
        if iFlag_verbose_in:
            logger.info(f'Feature processing summary:')
            logger.info(f'  - Total input features: {iFeature_index}')
            logger.info(f'  - Valid mesh cells created: {valid_mesh_cells}')
            logger.info(f'  - Invalid/skipped features: {invalid_geometry_count}')
            logger.info(f'  - Success rate: {((iFeature_index-invalid_geometry_count)/iFeature_index*100):.1f}%' if iFeature_index > 0 else '  - Success rate: 0%')
            # Report multipolygon handling statistics
            multipolygon_cells = valid_mesh_cells - (iFeature_index - invalid_geometry_count)
            if multipolygon_cells > 0:
                logger.info(f'  - Additional cells from multipolygons: {multipolygon_cells}')
                logger.info(f'  - Total mesh cells (including multipolygon parts): {valid_mesh_cells}')
        # Clean up dataset
        pDataset = None
        if not lons_list:
            logger.error('No valid polygon features found in mesh file')
            return None
        if iFlag_verbose_in:
            logger.info(f'Successfully processed {len(lons_list)} polygon features')
        # Calculate maximum vertices and pad coordinates efficiently
        try:
            if not lons_list:
                logger.error('No coordinate data found')
                return None
            max_vertices = max(len(coord) for coord in lons_list)
            if max_vertices == 0:
                logger.error('No vertices found in any polygon')
                return None
            nVertex_max = max_vertices
            if iFlag_verbose_in:
                logger.info(f'Maximum vertices per polygon: {max_vertices}')
            # Pre-allocate arrays for better memory efficiency
            num_polygons = len(lons_list)
            lons_padded = np.full((num_polygons, max_vertices), np.nan, dtype=np.float64)
            lats_padded = np.full((num_polygons, max_vertices), np.nan, dtype=np.float64)
            # Fill padded arrays efficiently
            for i, (lon_coords, lat_coords) in enumerate(zip(lons_list, lats_list)):
                # Ensure coordinates are numpy arrays with proper dtype
                lon_coords = np.asarray(lon_coords, dtype=np.float64)
                lat_coords = np.asarray(lat_coords, dtype=np.float64)
                # Validate coordinate data
                if len(lon_coords) != len(lat_coords):
                    logger.warning(f'Coordinate length mismatch in polygon {i}: lon={len(lon_coords)}, lat={len(lat_coords)}')
                    min_len = min(len(lon_coords), len(lat_coords))
                    lon_coords = lon_coords[:min_len]
                    lat_coords = lat_coords[:min_len]
                # Check for valid coordinate values
                if not (np.all(np.isfinite(lon_coords)) and np.all(np.isfinite(lat_coords))):
                    logger.warning(f'Invalid coordinates found in polygon {i}')
                    # Remove invalid coordinates
                    valid_mask = np.isfinite(lon_coords) & np.isfinite(lat_coords)
                    lon_coords = lon_coords[valid_mask]
                    lat_coords = lat_coords[valid_mask]
                coord_len = len(lon_coords)
                if coord_len > 0:
                    lons_padded[i, :coord_len] = lon_coords
                    lats_padded[i, :coord_len] = lat_coords
                else:
                    logger.warning(f'No valid coordinates remaining for polygon {i}')
            # Convert to the expected format for backward compatibility
            lons = lons_padded
            lats = lats_padded
        except Exception as e:
            logger.error(f'Error during coordinate padding: {str(e)}')
            logger.error(f'Traceback: {traceback.format_exc()}')
            return None
        # Calculate centroids efficiently using vectorized operations
        try:
            cell_lons_1d = []
            cell_lats_1d = []
            # Pre-allocate arrays for better performance
            cell_lons_1d = np.zeros(len(lons_list), dtype=np.float64)
            cell_lats_1d = np.zeros(len(lons_list), dtype=np.float64)
            for i in range(len(lons_list)):
                # Calculate centroid of each cell (ignoring NaN values)
                valid_mask = ~np.isnan(lons[i])
                if np.any(valid_mask):
                    valid_lons = lons[i][valid_mask]
                    valid_lats = lats[i][valid_mask]
                    # Use vectorized operations for better performance
                    centroid_lon = np.mean(valid_lons)
                    centroid_lat = np.mean(valid_lats)
                    # Validate centroid coordinates
                    if np.isfinite(centroid_lon) and np.isfinite(centroid_lat):
                        cell_lons_1d[i] = centroid_lon
                        cell_lats_1d[i] = centroid_lat
                    else:
                        logger.warning(f'Invalid centroid calculated for cell {i}: lon={centroid_lon}, lat={centroid_lat}')
                        # Use geometric center of bounding box as fallback
                        if len(valid_lons) > 0 and len(valid_lats) > 0:
                            cell_lons_1d[i] = (np.min(valid_lons) + np.max(valid_lons)) / 2.0
                            cell_lats_1d[i] = (np.min(valid_lats) + np.max(valid_lats)) / 2.0
                        else:
                            cell_lons_1d[i] = 0.0
                            cell_lats_1d[i] = 0.0
                else:
                    logger.warning(f'No valid coordinates found for cell {i}')
                    cell_lons_1d[i] = 0.0
                    cell_lats_1d[i] = 0.0
            if iFlag_verbose_in:
                logger.info(f'Calculated centroids for {len(cell_lons_1d)} cells')
            # Validate centroid ranges
            lon_range = (np.min(cell_lons_1d), np.max(cell_lons_1d))
            lat_range = (np.min(cell_lats_1d), np.max(cell_lats_1d))
            if not (-180 <= lon_range[0] <= 180 and -180 <= lon_range[1] <= 180):
                logger.warning(f'Longitude centroids outside valid range: {lon_range}')
            if not (-90 <= lat_range[0] <= 90 and -90 <= lat_range[1] <= 90):
                logger.warning(f'Latitude centroids outside valid range: {lat_range}')
        except Exception as e:
            logger.error(f'Error during centroid calculation: {str(e)}')
            return None
        # Extract unique vertices and connectivity
        try:
            if iFlag_verbose_in:
                logger.info('Extracting unique vertices and connectivity...')
            xv, yv, connectivity, vertex_to_index = extract_unique_vertices_and_connectivity(
                lons_list, lats_list
            )
            if xv is None or yv is None or connectivity is None:
                logger.error('Failed to extract unique vertices and connectivity')
                return None
            if iFlag_verbose_in:
                logger.info(f'Extracted {len(xv)} unique vertices')
                logger.info(f'Created connectivity matrix with shape: {connectivity.shape}')
        except Exception as e:
            logger.error(f'Error during vertex/connectivity extraction: {str(e)}')
            return None
        # Store results in class attributes
        aVertex_longititude = xv
        aVertex_latitude = yv
        aCenter_longititude = cell_lons_1d
        aCenter_latitude = cell_lats_1d
        aConnectivity = connectivity
        # Ensure aCellID matches the number of valid mesh cells
        if len(aCellID) != len(cell_lons_1d):
            logger.warning(f"aCellID length ({len(aCellID)}) doesn't match mesh cells ({len(cell_lons_1d)})")
            if len(aCellID) > len(cell_lons_1d):
                # Truncate aCellID to match mesh cells
                logger.warning("Truncating aCellID to match mesh cell count")
                aCellID = aCellID[:len(cell_lons_1d)]
            else:
                # Extend aCellID with sequential indices
                logger.warning("Extending aCellID with sequential indices to match mesh cell count")
                missing_count = len(cell_lons_1d) - len(aCellID)
                aCellID.extend(range(len(aCellID), len(aCellID) + missing_count))
        aCellID = np.array(aCellID)
        if iFlag_verbose_in:
            logger.info(f'Final aCellID array length: {len(aCellID)}')
        # Calculate and store area statistics
        if aArea_list:
            area_array = np.array(aArea_list)
            valid_areas = area_array[area_array > 0]  # Exclude zero areas from statistics
            if len(valid_areas) > 0:
                dArea_min = float(np.min(valid_areas))
                dArea_max = float(np.max(valid_areas))
                dArea_mean = float(np.mean(valid_areas))
                dArea_max = float(np.max(valid_areas))
                dArea_min = float(np.min(valid_areas))
                if iFlag_verbose_in:
                    logger.info(f'Mesh area statistics:')
                    logger.info(f'  - Min area: {dArea_min:.6f}')
                    logger.info(f'  - Max area: {dArea_max:.6f}')
                    logger.info(f'  - Mean area: {dArea_mean:.6f}')
            else:
                logger.warning('No valid polygon areas calculated')
                dArea_min = 0.0
                dArea_max = 0.0
                dArea_mean = 0.0
        # Enhanced validation of final results
        validation_passed = True
        if len(aVertex_longititude) == 0:
            logger.error('No unique vertices extracted')
            validation_passed = False
        if len(aCenter_longititude) != len(lons_list):
            logger.error(f'Centroid count mismatch: expected {len(lons_list)}, got {len(aCenter_longititude)}')
            validation_passed = False
        if aConnectivity is None or aConnectivity.size == 0:
            logger.error('Empty connectivity matrix')
            validation_passed = False
        # Validate connectivity indices
        if aConnectivity is not None:
            max_vertex_index = len(aVertex_longititude) - 1
            valid_connectivity = aConnectivity[aConnectivity >= 0]
            if len(valid_connectivity) > 0 and np.max(valid_connectivity) > max_vertex_index:
                logger.error('Connectivity matrix contains invalid vertex indices')
                validation_passed = False
        # Check for reasonable mesh bounds
        if len(aVertex_longititude) > 0:
            vertex_lon_range = (np.min(aVertex_longititude), np.max(aVertex_longititude))
            vertex_lat_range = (np.min(aVertex_latitude), np.max(aVertex_latitude))
            # Basic range reporting (detailed validation done by check_geometry_validity)
            if iFlag_verbose_in:
                logger.debug(f'Vertex longitude range: {vertex_lon_range}')
                logger.debug(f'Vertex latitude range: {vertex_lat_range}')
        if not validation_passed:
            logger.error('Mesh topology rebuild failed validation')
            return None
        if iFlag_verbose_in:
            logger.info('Mesh topology successfully rebuilt')
            logger.info(f'Final mesh statistics:')
            logger.info(f'  - Unique vertices: {len(aVertex_longititude)}')
            logger.info(f'  - Mesh cells: {len(aCenter_longititude)}')
            logger.info(f'  - Max vertices per cell: {nVertex_max}')
            logger.info(f'  - Connectivity shape: {aConnectivity.shape}')
            logger.info(f'  - Vertex longitude range: [{np.min(aVertex_longititude):.3f}, {np.max(aVertex_longititude):.3f}]')
            logger.info(f'  - Vertex latitude range: [{np.min(aVertex_latitude):.3f}, {np.max(aVertex_latitude):.3f}]')

        # Return comprehensive mesh topology information
        mesh_info = {
            'vertices_longitude': aVertex_longititude,
            'vertices_latitude': aVertex_latitude,
            'connectivity': aConnectivity,
            'cell_centroids_longitude': aCenter_longititude,
            'cell_centroids_latitude': aCenter_latitude,
            'cell_ids': aCellID,
            'area_min': dArea_min,
            'area_max': dArea_max,
            'area_mean': dArea_mean,
            'max_vertices_per_cell': nVertex_max,
            'num_cells': nFeatures,
            'num_polygns': len(aCenter_longititude),
            'num_vertices': len(aVertex_longititude),
            'success': True
        }

        return mesh_info
    except Exception as e:
        logger.error(f'Unexpected error in rebuild_mesh_topology: {str(e)}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return None

def normalize_longitude(lon):
    """
    Normalize longitude to [-180, 180] range by wrapping.
    Optimized using modular arithmetic instead of loops.

    Args:
        lon (float): Longitude value

    Returns:
        float: Normalized longitude in [-180, 180] range
    """
    # Use modular arithmetic for O(1) normalization
    normalized = ((lon + 180) % 360) - 180
    return normalized

def fix_coordinates(coords):
    """
    Recursively fix longitude values in coordinate arrays.
    Handles nested coordinate structures for different geometry types.

    Args:
        coords: Coordinate array (can be nested)

    Returns:
        Fixed coordinate array with normalized longitudes
    """
    if not coords:
        return coords

    # Check if this is a coordinate pair [lon, lat]
    if isinstance(coords[0], (int, float)):
        return [normalize_longitude(coords[0]), coords[1]]
    else:
        return [fix_coordinates(coord) for coord in coords]

def fix_longitude_range_gdal(geometry, in_place=False):
    """
    Fix longitude values using GDAL geometry operations.
    Normalizes longitude coordinates to [-180, 180] range.

    Args:
        geometry: OGR Geometry object
        in_place (bool): If True, modify geometry in place (faster)

    Returns:
        OGR Geometry object with normalized longitude coordinates
    """
    if geometry is None:
        return geometry

    # Optionally clone the geometry to avoid modifying the original
    fixed_geometry = geometry if in_place else geometry.Clone()

    # Get geometry type
    geom_type = fixed_geometry.GetGeometryName()

    if geom_type in ['POLYGON', 'MULTIPOLYGON', 'LINESTRING', 'MULTILINESTRING', 'POINT', 'MULTIPOINT']:
        # For complex geometries, iterate through all geometry parts
        if geom_type.startswith('MULTI') or geom_type == 'POLYGON':
            _fix_geometry_coordinates_recursive(fixed_geometry)
        else:
            # For simple geometries, fix coordinates directly using batch processing
            point_count = fixed_geometry.GetPointCount()
            if point_count > 0:
                # Process points in batches for better performance
                for i in range(point_count):
                    x, y, z = fixed_geometry.GetPoint(i)
                    # Normalize longitude using modular arithmetic
                    normalized_x = ((x + 180) % 360) - 180
                    # Avoid exact +180° by nudging to just under 180°
                    if abs(normalized_x - 180.0) < 1e-10:
                        normalized_x = 180.0 - 1e-8
                    fixed_geometry.SetPoint(i, normalized_x, y, z)

    return fixed_geometry

def _fix_geometry_coordinates_recursive(geometry):
    """
    Recursively fix coordinates in complex geometries.

    Args:
        geometry: OGR Geometry object to fix in-place
    """
    geom_count = geometry.GetGeometryCount()

    if geom_count > 0:
        # Recurse through sub-geometries
        for i in range(geom_count):
            sub_geom = geometry.GetGeometryRef(i)
            _fix_geometry_coordinates_recursive(sub_geom)
    else:
        # Fix coordinates in this geometry
        point_count = geometry.GetPointCount()
        for i in range(point_count):
            x, y, z = geometry.GetPoint(i)
            if abs(x - 180.0) < 1e-10:
                x = 180.0 - 1e-8  # Nudge to just under 180°
            if abs(x + 180.0) < 1e-10:
                x = -180.0 + 1e-8  # Nudge to just above -180°
            # Normalize longitude using modular arithmetic
            normalized_x = ((x + 180) % 360) - 180
            geometry.SetPoint(i, normalized_x, y, z)

def fix_mesh_longitude_range_and_idl_crossing(sFilename_in, sFilename_out, handle_idl_crossing=True):
    """
    Comprehensive GDAL-based function to fix longitude range issues and optionally handle
    International Date Line (IDL) crossing in vector files.

    This function combines the functionality of both fix_mesh_longitude_range_gdal and fix_idl_crossing
    into a single optimized function that can handle multiple layers and IDL crossing.

    Args:
        input_file (str): Path to input vector file (any GDAL-supported format)
        output_file (str): Path to output vector file
        handle_idl_crossing (bool): Whether to check and split polygons crossing the IDL (default: True)

    Returns:
        bool: True if successful, False otherwise
    """
    pDriver = None
    pDataset = None
    pDataset_out = None

    try:
        # Open source dataset
        pDriver = get_vector_driver_from_filename(sFilename_in)
        pDataset = pDriver.Open(sFilename_in, 0)
        if pDataset is None:
            logger.error(f"Could not open input file: {sFilename_in}")
            return False

        # Create output dataset
        pDriver_out = get_vector_driver_from_filename(sFilename_out)
        #delete output file if it already exists
        if os.path.exists(sFilename_out):
            pDriver_out.DeleteDataSource(sFilename_out)
        pDataset_out = pDriver_out.CreateDataSource(sFilename_out)
        if pDataset_out is None:
            logger.error(f"Could not create output file: {sFilename_out}")
            return False

        # Process each layer
        layer_count = pDataset.GetLayerCount()
        logger.info(f"Processing {layer_count} layer(s) from {sFilename_in}")

        total_processed = 0

        for layer_idx in range(layer_count):
            pLayer = pDataset.GetLayerByIndex(layer_idx)
            if pLayer is None:
                logger.warning(f"Could not get layer {layer_idx}")
                continue

            pLayerDefn = pLayer.GetLayerDefn()
            sSpatial_ref = pLayer.GetSpatialRef()
            layer_name = pLayer.GetName()
            nFeatures = pLayer.GetFeatureCount()

            logger.info(f"Processing layer '{layer_name}' with {nFeatures} features")

            # Create output layer with same schema
            pLayer_out = pDataset_out.CreateLayer(layer_name, sSpatial_ref, ogr.wkbUnknown)
            if pLayer_out is None:
                logger.error(f"Could not create output layer: {layer_name}")
                continue

            # Copy field definitions
            for iField in range(pLayerDefn.GetFieldCount()):
                field_defn = pLayerDefn.GetFieldDefn(iField)
                pLayer_out.CreateField(field_defn)

            # Process features with progress tracking
            processed_count = 0
            idl_crossing_count = 0

            for pFeature in pLayer:
                try:
                    geometry = pFeature.GetGeometryRef()
                    if geometry is None:
                        logger.warning(f"Feature ID {pFeature.GetFID()} has no geometry, skipping...")
                        continue

                    geometry_type = geometry.GetGeometryName()

                    #check whether geometry contains poles
                    aCoord_origin = get_geometry_coordinates(geometry)
                    if np.min(np.abs(aCoord_origin[:, 1])) >88.0:
                        continue
                    # Fix longitude coordinates using GDAL
                    fixed_geometry = fix_longitude_range_gdal(geometry)
                    # Handle IDL crossing for polygon geometries if requested
                    if handle_idl_crossing and geometry_type in ['POLYGON', 'MULTIPOLYGON']:
                        # Check for IDL crossing after longitude normalization
                        aCoord = get_geometry_coordinates(fixed_geometry)
                        bCross_idl, aCoord_updated = check_cross_international_date_line_polygon(aCoord)

                        if bCross_idl:
                            idl_crossing_count += 1
                            logger.info(f"Feature ID {pFeature.GetFID()} crosses the International Date Line. Splitting...")

                            [eastern_polygon, western_polygon] = split_international_date_line_polygon_coordinates(aCoord)

                            # Create a multipolygon geometry
                            pGeometry_multi = ogr.Geometry(ogr.wkbMultiPolygon)

                            # Create eastern polygon
                            pPolygon_eastern = ogr.Geometry(ogr.wkbPolygon)
                            pLinearRing_eastern = ogr.Geometry(ogr.wkbLinearRing)
                            for coord in eastern_polygon:
                                pLinearRing_eastern.AddPoint(coord[0], coord[1])
                            pLinearRing_eastern.CloseRings()
                            pPolygon_eastern.AddGeometry(pLinearRing_eastern)
                            pGeometry_multi.AddGeometry(pPolygon_eastern)

                            # Create western polygon
                            pPolygon_western = ogr.Geometry(ogr.wkbPolygon)
                            pLinearRing_western = ogr.Geometry(ogr.wkbLinearRing)
                            for coord in western_polygon:
                                pLinearRing_western.AddPoint(coord[0], coord[1])
                            pLinearRing_western.CloseRings()
                            pPolygon_western.AddGeometry(pLinearRing_western)
                            pGeometry_multi.AddGeometry(pPolygon_western)

                            final_geometry = pGeometry_multi
                        else:
                            if aCoord_updated is not None:
                                # Update fixed_geometry with adjusted coordinates to create a polygon, usually because of IDL
                                fixed_geometry = create_geometry_from_coordinates(aCoord_updated, geometry_type)
                            else:
                                final_geometry = fixed_geometry
                    else:
                        final_geometry = fixed_geometry

                    # Create output feature
                    pFeature_out = ogr.Feature(pLayer_out.GetLayerDefn())
                    pFeature_out.SetGeometry(final_geometry)

                    # Copy all field values
                    for iField in range(pLayerDefn.GetFieldCount()):
                        sField_name = pLayerDefn.GetFieldDefn(iField).GetName()
                        pFeature_out.SetField(sField_name, pFeature.GetField(sField_name))

                    pLayer_out.CreateFeature(pFeature_out)
                    pFeature_out = None

                    processed_count += 1
                    if processed_count % 1000 == 0:
                        logger.info(f"Processed {processed_count}/{nFeatures} features in layer '{layer_name}'...")

                except Exception as e:
                    logger.error(f"Error processing feature ID {pFeature.GetFID()} in layer '{layer_name}': {str(e)}")
                    continue

            if handle_idl_crossing and idl_crossing_count > 0:
                logger.info(f"Layer '{layer_name}': {processed_count} features processed, {idl_crossing_count} IDL crossings handled")
            else:
                logger.info(f"Layer '{layer_name}': {processed_count} features processed")

            total_processed += processed_count

        # Cleanup
        pDataset_out.FlushCache()
        pDataset_out = None
        pDataset = None

        logger.info(f"Successfully processed {total_processed} total features and created fixed file: {sFilename_out}")
        return True

    except Exception as e:
        logger.error(f"Error in fix_mesh_longitude_range_and_idl: {str(e)}")
        return False

    finally:
        # Ensure proper cleanup
        if pDataset_out is not None:
            pDataset_out.FlushCache()
            pDataset_out = None
        if pDataset is not None:
            pDataset = None

def create_geometry_from_coordinates(aCoord, geometry_type):
    """
    Create an OGR Geometry object from coordinate array based on specified geometry type.

    Args:
        aCoord (np.ndarray): Array of coordinates
        geometry_type (str): Type of geometry ('POLYGON', 'LINESTRING', etc.)

    Returns:
        OGR Geometry object
    """
    if geometry_type == 'POLYGON':
        pPolygon = ogr.Geometry(ogr.wkbPolygon)
        pLinearRing = ogr.Geometry(ogr.wkbLinearRing)
        for coord in aCoord:
            pLinearRing.AddPoint(coord[0], coord[1])
        pLinearRing.CloseRings()
        pPolygon.AddGeometry(pLinearRing)
        return pPolygon
    elif geometry_type == 'LINESTRING':
        pLineString = ogr.Geometry(ogr.wkbLineString)
        for coord in aCoord:
            pLineString.AddPoint(coord[0], coord[1])
        return pLineString
    elif geometry_type == 'POINT':
        pPoint = ogr.Geometry(ogr.wkbPoint)
        pPoint.AddPoint(aCoord[0][0], aCoord[0][1])
        return pPoint
    else:
        logger.error(f"Unsupported geometry type for creation: {geometry_type}")
        return None

