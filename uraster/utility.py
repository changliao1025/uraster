
import os
import logging
import traceback
from typing import Optional, Tuple, List, Dict, Any, Union
from osgeo import gdal, ogr

from uraster.classes.sraster import sraster
gdal.UseExceptions()
from pyearth.gis.geometry.international_date_line_utility import split_international_date_line_polygon_coordinates, check_cross_international_date_line_polygon
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
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

def _validate_geometry(geometry: Any, feature_id: Union[int, str], attempt_fix: bool = True) -> bool:
    """
    Validate and optionally fix OGR geometry objects.

    Args:
        geometry: OGR geometry object to validate
        feature_id: Feature identifier for logging
        attempt_fix: Whether to attempt fixing invalid geometries

    Returns:
        bool: True if geometry is valid (or was successfully fixed), False otherwise
    """
    if geometry is None:
        logger.warning(f"Geometry is None for feature {feature_id}")
        return False

    if geometry.IsValid():
        return True

    if not attempt_fix:
        logger.warning(f"Invalid geometry for feature {feature_id}")
        return False

    try:
        # Attempt to fix invalid geometry using buffer(0) technique
        logger.warning(
            f"Invalid geometry for feature {feature_id}, attempting to fix")
        fixed_geometry = geometry.Buffer(0)

        if fixed_geometry and fixed_geometry.IsValid():
            # Copy the fixed geometry back to the original
            geometry.ImportFromWkt(fixed_geometry.ExportToWkt())
            logger.debug(
                f"Successfully fixed geometry for feature {feature_id}")
            return True
        else:
            logger.warning(
                f"Cannot fix invalid geometry for feature {feature_id}")
            return False

    except Exception as e:
        logger.warning(
            f"Error attempting to fix geometry for feature {feature_id}: {e}")
        return False

def _log_memory_usage(stage: str, iFlag_verbose: bool = False) -> None:
    """
    Log current memory usage if psutil is available.

    Args:
        stage: Description of the current processing stage
        iFlag_verbose: Whether to log memory information
    """
    if not PSUTIL_AVAILABLE or not iFlag_verbose:
        return

    try:
        process = psutil.Process()
        memory_info = process.memory_info()
        memory_mb = memory_info.rss / (1024 * 1024)
        logger.info(f"Memory usage at {stage}: {memory_mb:.1f} MB")
    except Exception as e:
        logger.debug(f"Could not get memory usage: {e}")

def get_polygon_list(
    sFilename_source_mesh: str,
    iFlag_verbose_in: bool = False
) -> Optional[Tuple[List[Tuple[Union[int, str], str]], List[float], Optional[str]]]:
    """
    Extract polygon geometries and areas from mesh vector file.

    Processes mesh features, handles International Date Line (IDL) crossing polygons,
    and returns polygon WKT strings with associated areas and projection information.

    Args:
        sFilename_source_mesh (str): Path to the source mesh vector file
        iFlag_verbose (bool, optional): If True, print detailed progress messages.
            Default is False.

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
                # Get geometry and validate
                pPolygon = pFeature_mesh.GetGeometryRef()
                if pPolygon is None:
                    logger.warning(f"Feature {i} has no geometry, skipping")
                    error_count += 1
                    i += 1
                    continue

                if not pPolygon.IsValid():
                    logger.warning(
                        f"Feature {i} has invalid geometry, attempting to fix")
                    pPolygon = pPolygon.Buffer(0)  # Attempt to fix
                    if not pPolygon.IsValid():
                        logger.warning(
                            f"Cannot fix invalid geometry for feature {i}, skipping")
                        error_count += 1
                        i += 1
                        continue

                sGeometry_type = pPolygon.GetGeometryName()

                # Read cellid from current feature with error handling
                try:
                    current_cellid = pFeature_mesh.GetField('cellid')
                    if current_cellid is None:
                        current_cellid = i  # Use feature index as fallback
                except Exception as field_error:
                    logger.warning(
                        f"Error reading cellid for feature {i}: {field_error}")
                    current_cellid = i

                if sGeometry_type == "POLYGON":
                    try:
                        aCoord = get_geometry_coordinates(pPolygon)
                        if aCoord is None or len(aCoord) < 3:
                            logger.warning(
                                f"Invalid coordinates for polygon feature {i}, skipping")
                            error_count += 1
                            i += 1
                            continue

                        # Check whether geometry crosses the International Date Line
                        if check_cross_international_date_line_polygon(aCoord):
                            dArea = 0.0
                            if iFlag_verbose_in:
                                logger.info(
                                    f'Feature {i} crosses the international date line, splitting into multiple parts.')

                            # Create multipolygon to handle IDL crossing
                            pMultipolygon = ogr.Geometry(ogr.wkbMultiPolygon)
                            aCoord_gcs_split = split_international_date_line_polygon_coordinates(
                                aCoord)

                            for aCoord_gcs in aCoord_gcs_split:
                                try:
                                    dArea += calculate_polygon_area(
                                        aCoord_gcs[:, 0], aCoord_gcs[:, 1])

                                    # Create polygon from coordinates
                                    ring = ogr.Geometry(ogr.wkbLinearRing)
                                    for iCoord in range(aCoord_gcs.shape[0]):
                                        ring.AddPoint(
                                            aCoord_gcs[iCoord, 0], aCoord_gcs[iCoord, 1])

                                    ring.CloseRings()
                                    polygon_part = ogr.Geometry(ogr.wkbPolygon)
                                    polygon_part.AddGeometry(ring)

                                    # Validate polygon part
                                    if not polygon_part.IsValid():
                                        polygon_part = polygon_part.Buffer(
                                            0)  # Attempt to fix
                                        if not polygon_part.IsValid():
                                            logger.warning(
                                                f'Polygon part of feature {i} is invalid after IDL splitting')
                                            continue

                                    pMultipolygon.AddGeometry(polygon_part)

                                except Exception as part_error:
                                    logger.warning(
                                        f"Error processing IDL part for feature {i}: {part_error}")
                                    continue

                            if pMultipolygon.GetGeometryCount() > 0:
                                wkt = pMultipolygon.ExportToWkt()
                                aPolygon.append((current_cellid, wkt))
                                aArea.append(dArea)
                                processed_count += 1
                            else:
                                logger.warning(
                                    f"No valid parts created for IDL feature {i}")
                                error_count += 1
                        else:
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