# Define a class named 'uraster'
import os
import logging
import traceback
import numpy as np
from osgeo import gdal, ogr, osr
gdal.UseExceptions()
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
from pyearth.gis.geometry.extract_unique_vertices_and_connectivity import extract_unique_vertices_and_connectivity
from uraster.classes.sraster import sraster
from uraster.classes import _visual
from uraster.operation import extract

# Set up logging for crash detection
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)
crs = "EPSG:4326"
pDriver_geojson = ogr.GetDriverByName('GeoJSON')
pDriver_shp = ogr.GetDriverByName('ESRI Shapefile')

class uraster:
    """
    Unstructured raster processing class for zonal statistics on mesh geometries.

    Handles complex scenarios including:
    - International Date Line (IDL) crossing polygons
    - Parallel processing for large datasets
    - Multiple raster formats and coordinate systems
    - Comprehensive error handling and crash detection
    """

    def __init__(self, aConfig=None):
        """
        Initialize uraster instance.

        Args:
            aConfig (dict, optional): Configuration dictionary with keys:
                - iFlag_remap_method (int): Remap method (1=nearest, 2=nearest, 3=weighted average)
                - sFilename_source_mesh (str): Source mesh file path
                - sFilename_target_mesh (str): Target mesh file path
                - aFilename_source_raster (list): List of source raster file paths
        """
        # Default configuration
        if aConfig is None:
            aConfig = {}

        # Processing flags and resolutions
        self.iFlag_global = None
        self.iFlag_remap_method = aConfig.get('iFlag_remap_method', 1)  # Default to nearest neighbor
        self.dResolution_raster = None
        self.dResolution_uraster = None

        # File paths
        self.sFilename_source_mesh = aConfig.get('sFilename_source_mesh', None)
        self.sField_unique_id = aConfig.get('sField_unique_id', None)
        self.sFilename_target_mesh = aConfig.get('sFilename_target_mesh', None)
        self.aFilename_source_raster = aConfig.get('aFilename_source_raster', [])

        # Cell counts
        self.nCell = -1
        self.nCell_source = -1
        self.nCell_target = -1
        self.nVertex_max = 0  # Will be calculated dynamically

        # Mesh topology data
        self.aVertex_longititude = None
        self.aVertex_latitude = None
        self.aCenter_longititude = None
        self.aCenter_latitude = None
        self.aConnectivity = None
        self.aCellID = None

        # Mesh area statistics
        self.dArea_min = None
        self.dArea_max = None
        self.dArea_mean = None

        # Resolution comparison threshold (ratio of mesh to raster resolution)
        # If mesh cells are within this factor of raster resolution, use weighted averaging
        self.dResolution_ratio_threshold = 3.0  # mesh resolution < 3x raster resolution triggers weighted avg

        # Validate configuration
        if self.iFlag_remap_method not in [1, 2, 3]:
            logger.warning(f"Invalid remap method {self.iFlag_remap_method}, defaulting to 1 (nearest neighbor)")
            self.iFlag_remap_method = 1

    def setup(self, iFlag_verbose=False):
        """
        Initialize and validate the uraster configuration.
        Checks raster files and mesh file for existence and validity.

        Args:
            iFlag_verbose (bool, optional): If True, print detailed progress messages.
                If False, only print error messages. Default is False.

        Returns:
            bool: True if setup successful, False otherwise
        """
        raster_check = self.check_raster_files(iFlag_verbose=iFlag_verbose)
        mesh_check = self.check_mesh_file(iFlag_verbose=iFlag_verbose)

        return raster_check is not None and mesh_check is not None

    def check_raster_files(self, aFilename_source_raster_in=None, iFlag_verbose=False):
        """
        Validate and prepare input raster files, converting to WGS84 if needed.

        Performs comprehensive validation of raster files including:
        - File existence and readability
        - Valid GDAL raster format
        - Coordinate system compatibility
        - Data integrity checks

        Args:
            aFilename_source_raster_in (list, optional): List of raster file paths.
                If None, uses self.aFilename_source_raster
            iFlag_verbose (bool, optional): If True, print detailed progress messages.
                If False, only print error messages. Default is False.

        Returns:
            list: List of WGS84 raster file paths, or None if validation fails

        Note:
            - Non-WGS84 rasters are automatically converted and cached
            - All rasters must be valid and readable for processing to continue
        """
        # Determine input raster list
        if aFilename_source_raster_in is None:
            aFilename_source_raster = self.aFilename_source_raster
        else:
            aFilename_source_raster = aFilename_source_raster_in

        # Validate input list
        if not aFilename_source_raster:
            logger.error('No raster files provided for validation')
            return None

        if not isinstance(aFilename_source_raster, (list, tuple)):
            logger.error(f'Raster files must be provided as a list, got {type(aFilename_source_raster).__name__}')
            return None

        if iFlag_verbose:
            logger.info(f'Validating {len(aFilename_source_raster)} raster file(s)...')

        # Phase 1: Check file existence and readability
        for idx, sFilename_raster_in in enumerate(aFilename_source_raster, 1):
            if not isinstance(sFilename_raster_in, str):
                logger.error(f'Raster file path must be a string, got {type(sFilename_raster_in).__name__} at index {idx}')
                return None

            if not sFilename_raster_in.strip():
                logger.error(f'Empty raster file path at index {idx}')
                return None

            if not os.path.exists(sFilename_raster_in):
                logger.error(f'Raster file does not exist: {sFilename_raster_in}')
                return None

            if not os.path.isfile(sFilename_raster_in):
                logger.error(f'Path is not a file: {sFilename_raster_in}')
                return None

            # Check file permissions
            if not os.access(sFilename_raster_in, os.R_OK):
                logger.error(f'Raster file is not readable: {sFilename_raster_in}')
                return None

            # Quick GDAL format validation
            try:
                pDataset_test = gdal.Open(sFilename_raster_in, gdal.GA_ReadOnly)
                if pDataset_test is None:
                    logger.error(f'GDAL cannot open raster file: {sFilename_raster_in}')
                    return None
                pDataset_test = None  # Close dataset
            except Exception as e:
                logger.error(f'Error opening raster with GDAL: {sFilename_raster_in}: {e}')
                return None

        if iFlag_verbose:
            logger.info('All raster files exist and are readable')

        # Phase 2: Process and convert rasters to WGS84
        aFilename_source_raster_out = []

        # Create WGS84 spatial reference for comparison
        pSpatialRef_wgs84 = None
        try:
            pSpatialRef_wgs84 = osr.SpatialReference()
            pSpatialRef_wgs84.ImportFromEPSG(4326)
            wkt_wgs84 = pSpatialRef_wgs84.ExportToWkt()
        except Exception as e:
            logger.error(f'Failed to create WGS84 spatial reference: {e}')
            return None
        finally:
            # Clean up spatial reference object
            if pSpatialRef_wgs84 is not None:
                pSpatialRef_wgs84 = None

        # Process each raster file
        for idx, sFilename_raster_in in enumerate(aFilename_source_raster, 1):
            if iFlag_verbose:
                logger.info(f'Processing raster {idx}/{len(aFilename_source_raster)}: {os.path.basename(sFilename_raster_in)}')

            try:
                # Create sraster instance and read metadata
                pRaster = sraster(sFilename_raster_in)
                pRaster.read_metadata()

                # Validate critical metadata
                if pRaster.pSpatialRef_wkt is None:
                    logger.error(f'Raster has no spatial reference: {sFilename_raster_in}')
                    return None

                if pRaster.nrow is None or pRaster.ncolumn is None:
                    logger.error(f'Invalid raster dimensions: {sFilename_raster_in}')
                    return None

                if pRaster.nrow <= 0 or pRaster.ncolumn <= 0:
                    logger.error(f'Raster has invalid dimensions ({pRaster.nrow}x{pRaster.ncolumn}): {sFilename_raster_in}')
                    return None

                # Check if coordinate system matches WGS84
                if pRaster.pSpatialRef_wkt == wkt_wgs84:
                    if iFlag_verbose:
                        logger.info(f'  ✓ Already in WGS84 (EPSG:4326)')
                    aFilename_source_raster_out.append(sFilename_raster_in)
                else:
                    # Convert to WGS84
                    if iFlag_verbose:
                        logger.info(f'  → Converting to WGS84 from {pRaster.pSpatialRef.GetName() if pRaster.pSpatialRef else "unknown CRS"}')
                    try:
                        pRaster_wgs84 = pRaster.convert_to_wgs84()

                        if pRaster_wgs84 is None or not hasattr(pRaster_wgs84, 'sFilename'):
                            logger.error(f'Conversion to WGS84 failed: {sFilename_raster_in}')
                            return None

                        if not os.path.exists(pRaster_wgs84.sFilename):
                            logger.error(f'Converted WGS84 file not found: {pRaster_wgs84.sFilename}')
                            return None

                        if iFlag_verbose:
                            logger.info(f'  ✓ Converted to: {pRaster_wgs84.sFilename}')
                        aFilename_source_raster_out.append(pRaster_wgs84.sFilename)

                    except Exception as e:
                        logger.error(f'Error during WGS84 conversion: {sFilename_raster_in}: {e}')
                        logger.error(f'Traceback: {traceback.format_exc()}')
                        return None

                # Log raster summary
                if iFlag_verbose:
                    logger.debug(f'  - Dimensions: {pRaster.nrow} x {pRaster.ncolumn} pixels')
                    logger.debug(f'  - Data type: {pRaster.eType}')
                    if hasattr(pRaster, 'dNoData'):
                        logger.debug(f'  - NoData value: {pRaster.dNoData}')

                pRaster.pSpatialRef = None  # Clean up spatial reference

            except AttributeError as e:
                logger.error(f'Missing expected attribute in sraster: {sFilename_raster_in}: {e}')
                logger.error(f'Ensure sraster class has all required methods and attributes')
                return None

            except Exception as e:
                logger.error(f'Unexpected error processing raster {sFilename_raster_in}: {e}')
                logger.error(f'Error type: {type(e).__name__}')
                logger.error(f'Traceback: {traceback.format_exc()}')
                return None

        # Final validation
        if len(aFilename_source_raster_out) != len(aFilename_source_raster):
            logger.error(f'Output count mismatch: expected {len(aFilename_source_raster)}, got {len(aFilename_source_raster_out)}')
            return None

        if iFlag_verbose:
            logger.info(f'Successfully validated and prepared {len(aFilename_source_raster_out)} raster file(s)')

        return aFilename_source_raster_out

    def check_mesh_file(self, iFlag_verbose=False):
        """
        Check if the source mesh file exists and build its topology.

        Args:
            iFlag_verbose (bool, optional): If True, print detailed progress messages.
                If False, only print error messages. Default is False.

        Returns:
            tuple or None: (vertices_lon, vertices_lat, connectivity) if successful, None otherwise
        """
        if not self.sFilename_source_mesh:
            logger.error("No source mesh filename provided")
            return None

        if not os.path.exists(self.sFilename_source_mesh):
            logger.error(f"Source mesh file does not exist: {self.sFilename_source_mesh}")
            return None

        return self.rebuild_mesh_topology(iFlag_verbose=iFlag_verbose)

    def _get_geometry_type_name(self, geometry_type):
        """
        Convert OGR geometry type integer to readable string name.

        Handles both standard and 3D/Z-flagged geometry types.

        Args:
            geometry_type (int): OGR geometry type constant

        Returns:
            str: Human-readable geometry type name (e.g., "wkbPolygon")
        """
        geometry_types = {
            ogr.wkbUnknown: "wkbUnknown",
            ogr.wkbPoint: "wkbPoint",
            ogr.wkbLineString: "wkbLineString",
            ogr.wkbPolygon: "wkbPolygon",
            ogr.wkbMultiPoint: "wkbMultiPoint",
            ogr.wkbMultiLineString: "wkbMultiLineString",
            ogr.wkbMultiPolygon: "wkbMultiPolygon",
            ogr.wkbGeometryCollection: "wkbGeometryCollection"
        }

        # Direct match
        if geometry_type in geometry_types:
            return geometry_types[geometry_type]

        # Check base type (removes 3D/Z flags)
        base_type = geometry_type & 0xFF
        for const_val, name in geometry_types.items():
            if (const_val & 0xFF) == base_type:
                return f"{name} (with flags)"

        return f"Unknown geometry type: {geometry_type}"

    def rebuild_mesh_topology(self, iFlag_verbose=False):
        """
        Rebuild mesh topology from source mesh file by extracting vertices,
        connectivity, and centroids for unstructured mesh processing.

        Args:
            iFlag_verbose (bool, optional): If True, print detailed progress messages.
                If False, only print error messages. Default is False.

        Returns:
            tuple: (vertices_longitude, vertices_latitude, connectivity) or None on failure
        """

        try:
            # Open the input data source
            pDataset = ogr.Open(self.sFilename_source_mesh, 0)  # Read-only
            if pDataset is None:
                logger.error(f'Failed to open file: {self.sFilename_source_mesh}')
                return None

            if iFlag_verbose:
                logger.info(f'Successfully opened mesh file: {self.sFilename_source_mesh}')

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
            self.nCell_source = nFeatures #if there is no invalid features, this is the number of cells
            iFieldCount = pLayerDefn.GetFieldCount()

            if nFeatures == 0:
                logger.warning('Layer contains no features.')
                pDataset = None
                return None

            aCellID= []  # Will be populated dynamically as features are processed

            if iFlag_verbose:
                logger.info(f'Processing {nFeatures} features with {iFieldCount} fields')

            # Get the first field name (assuming it contains the data variable)
            if self.sField_unique_id is None:
                sVariable = pLayerDefn.GetFieldDefn(0).GetName() if iFieldCount > 0 else None
                #search whether it has a field name used for id
                #for idx, pFeature_base in enumerate(pLayer):
                #    fid = pFeature_base.GetFID()
                #sVariable = None
            else:
                sVariable = self.sField_unique_id

            # Initialize lists for storing geometry data
            lons_list = []
            lats_list = []
            data_list = []
            area_list = []

            # Process features with enhanced error handling
            pLayer.ResetReading()
            iFeature_index = 0
            invalid_geometry_count = 0

            for pFeature in pLayer:
                if pFeature is None:
                    continue

                pGeometry = pFeature.GetGeometryRef()
                if pGeometry is None:
                    logger.warning(f'Feature {iFeature_index} has no geometry, skipping')
                    iFeature_index += 1
                    invalid_geometry_count += 1
                    continue

                sGeometry_type = pGeometry.GetGeometryName()
                if sGeometry_type == 'POLYGON':
                    try:
                        # Validate geometry before processing
                        if not pGeometry.IsValid():
                            logger.warning(f'Feature {iFeature_index} has invalid geometry')
                            invalid_geometry_count += 1
                            print(pGeometry.ExportToWkt())
                            continue
                        # Get coordinates of the polygon
                        aCoord = get_geometry_coordinates(pGeometry)
                        if aCoord is not None and len(aCoord) > 0:
                            # Validate coordinate bounds
                            lons = aCoord[:, 0]
                            lats = aCoord[:, 1]
                            # Check for reasonable coordinate ranges
                            if (np.any(lons < -180) or np.any(lons > 180) or
                                np.any(lats < -90) or np.any(lats > 90)):
                                logger.warning(f'Feature {iFeature_index} has coordinates outside valid range')

                            # Check for minimum polygon area (avoid degenerate polygons)
                            if len(aCoord) < 3:
                                logger.warning(f'Feature {iFeature_index} has fewer than 3 vertices, skipping')
                                iFeature_index += 1
                                invalid_geometry_count += 1
                                continue

                            lons_list.append(lons)
                            lats_list.append(lats)

                            # Calculate polygon area
                            try:
                                dArea = calculate_polygon_area(lons, lats)
                                area_list.append(dArea)
                            except Exception as area_error:
                                logger.warning(f'Could not calculate area for feature {iFeature_index}: {area_error}')
                                area_list.append(0.0)

                            # Get field data if available
                            if sVariable:
                                try:
                                    field_value = pFeature.GetField(sVariable)
                                    # Handle different field types
                                    if field_value is not None:
                                        data_list.append(int(field_value))
                                    else:
                                        data_list.append(0)

                                    aCellID.append(int(field_value) if field_value is not None else iFeature_index)
                                except (ValueError, TypeError) as e:
                                    logger.warning(f'Could not convert field value for feature {iFeature_index}: {e}')
                                    data_list.append(iFeature_index)
                                    aCellID.append(iFeature_index)
                            else:
                                data_list.append(iFeature_index)  # Use feature index as default
                                aCellID.append(iFeature_index)
                        else:
                            logger.warning(f'Failed to extract coordinates from feature {iFeature_index}')
                            invalid_geometry_count += 1

                    except Exception as e:
                        logger.warning(f'Error processing feature {iFeature_index}: {str(e)}')
                        invalid_geometry_count += 1

                elif sGeometry_type == 'MULTIPOLYGON':
                    try:
                        # Process multipolygon by extracting each constituent polygon
                        if iFlag_verbose:
                            logger.info(f'Processing multipolygon feature {iFeature_index} with {pGeometry.GetGeometryCount()} parts')

                        multipolygon_processed = False

                        for iPart in range(pGeometry.GetGeometryCount()):
                            pPolygon_part = pGeometry.GetGeometryRef(iPart)

                            if pPolygon_part is None:
                                logger.warning(f'Multipolygon part {iPart} is None in feature {iFeature_index}')
                                continue

                            if not pPolygon_part.IsValid():
                                logger.warning(f'Multipolygon part {iPart} has invalid geometry in feature {iFeature_index}')
                                continue

                            # Get coordinates of the polygon part
                            aCoord_part = get_geometry_coordinates(pPolygon_part)
                            if aCoord_part is not None and len(aCoord_part) > 0:
                                # Validate coordinate bounds for this part
                                lons_part = aCoord_part[:, 0]
                                lats_part = aCoord_part[:, 1]

                                # Check for reasonable coordinate ranges
                                if (np.any(lons_part < -180) or np.any(lons_part > 180) or
                                    np.any(lats_part < -90) or np.any(lats_part > 90)):
                                    logger.warning(f'Multipolygon part {iPart} in feature {iFeature_index} has coordinates outside valid range')

                                # Check for minimum polygon area (avoid degenerate polygons)
                                if len(aCoord_part) < 3:
                                    logger.warning(f'Multipolygon part {iPart} in feature {iFeature_index} has fewer than 3 vertices, skipping part')
                                    continue

                                lons_list.append(lons_part)
                                lats_list.append(lats_part)

                                # Calculate polygon area for this part
                                try:
                                    dArea_part = calculate_polygon_area(lons_part, lats_part)
                                    area_list.append(dArea_part)
                                except Exception as area_error:
                                    logger.warning(f'Could not calculate area for multipolygon part {iPart} in feature {iFeature_index}: {area_error}')
                                    area_list.append(0.0)

                                # For multipolygon, use the original feature index for all parts
                                # but track that this is a multipolygon part
                                if sVariable:
                                    try:
                                        field_value = pFeature.GetField(sVariable)
                                        if field_value is not None:
                                            data_list.append(int(field_value))
                                        else:
                                            data_list.append(0)
                                        # Maintain original aCellID for multipolygon features
                                        # Each part gets the same CellID as the original feature
                                        aCellID.append(int(field_value) if field_value is not None else iFeature_index)
                                    except (ValueError, TypeError) as e:
                                        logger.warning(f'Could not convert field value for multipolygon part {iPart} in feature {iFeature_index}: {e}')
                                        data_list.append(iFeature_index)
                                        aCellID.append(iFeature_index)
                                else:
                                    data_list.append(iFeature_index)
                                    aCellID.append(iFeature_index)

                                multipolygon_processed = True
                            else:
                                logger.warning(f'Failed to extract coordinates from multipolygon part {iPart} in feature {iFeature_index}')

                        if not multipolygon_processed:
                            logger.warning(f'No valid parts found in multipolygon feature {iFeature_index}')
                            invalid_geometry_count += 1
                        else:
                            if iFlag_verbose:
                                logger.info(f'Successfully processed multipolygon feature {iFeature_index}')

                    except Exception as e:
                        logger.warning(f'Error processing multipolygon feature {iFeature_index}: {str(e)}')
                        invalid_geometry_count += 1

                elif sGeometry_type in ['POINT', 'LINESTRING']:
                    logger.warning(f'Geometry type {sGeometry_type} not supported in feature {iFeature_index}, skipping')
                    invalid_geometry_count += 1
                else:
                    logger.warning(f'Unknown geometry type {sGeometry_type} in feature {iFeature_index}, skipping')
                    invalid_geometry_count += 1

                iFeature_index += 1

            # Report processing statistics
            valid_mesh_cells = len(lons_list)
            if iFlag_verbose:
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

            if iFlag_verbose:
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

                self.nVertex_max = max_vertices
                if iFlag_verbose:
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

                if iFlag_verbose:
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
                if iFlag_verbose:
                    logger.info('Extracting unique vertices and connectivity...')
                xv, yv, connectivity, vertex_to_index = extract_unique_vertices_and_connectivity(
                    lons_list, lats_list
                )

                if xv is None or yv is None or connectivity is None:
                    logger.error('Failed to extract unique vertices and connectivity')
                    return None

                if iFlag_verbose:
                    logger.info(f'Extracted {len(xv)} unique vertices')
                    logger.info(f'Created connectivity matrix with shape: {connectivity.shape}')

            except Exception as e:
                logger.error(f'Error during vertex/connectivity extraction: {str(e)}')
                return None

            # Store results in class attributes
            self.aVertex_longititude = xv
            self.aVertex_latitude = yv
            self.aCenter_longititude = cell_lons_1d
            self.aCenter_latitude = cell_lats_1d
            self.aConnectivity = connectivity

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

            self.aCellID = np.array(aCellID)

            if iFlag_verbose:
                logger.info(f'Final aCellID array length: {len(self.aCellID)}')
                logger.info(f'aCellID range: [{np.min(self.aCellID)}, {np.max(self.aCellID)}]')
            # Calculate and store area statistics
            if area_list:
                area_array = np.array(area_list)
                valid_areas = area_array[area_array > 0]  # Exclude zero areas from statistics
                if len(valid_areas) > 0:
                    self.dArea_min = float(np.min(valid_areas))
                    self.dArea_max = float(np.max(valid_areas))
                    self.dArea_mean = float(np.mean(valid_areas))
                    if iFlag_verbose:
                        logger.info(f'Mesh area statistics:')
                        logger.info(f'  - Min area: {self.dArea_min:.6f}')
                        logger.info(f'  - Max area: {self.dArea_max:.6f}')
                        logger.info(f'  - Mean area: {self.dArea_mean:.6f}')
                else:
                    logger.warning('No valid polygon areas calculated')
                    self.dArea_min = 0.0
                    self.dArea_max = 0.0
                    self.dArea_mean = 0.0

            # Enhanced validation of final results
            validation_passed = True

            if len(self.aVertex_longititude) == 0:
                logger.error('No unique vertices extracted')
                validation_passed = False

            if len(self.aCenter_longititude) != len(lons_list):
                logger.error(f'Centroid count mismatch: expected {len(lons_list)}, got {len(self.aCenter_longititude)}')
                validation_passed = False

            if self.aConnectivity is None or self.aConnectivity.size == 0:
                logger.error('Empty connectivity matrix')
                validation_passed = False

            # Validate connectivity indices
            if self.aConnectivity is not None:
                max_vertex_index = len(self.aVertex_longititude) - 1
                valid_connectivity = self.aConnectivity[self.aConnectivity >= 0]
                if len(valid_connectivity) > 0 and np.max(valid_connectivity) > max_vertex_index:
                    logger.error('Connectivity matrix contains invalid vertex indices')
                    validation_passed = False

            # Check for reasonable mesh bounds
            if len(self.aVertex_longititude) > 0:
                vertex_lon_range = (np.min(self.aVertex_longititude), np.max(self.aVertex_longititude))
                vertex_lat_range = (np.min(self.aVertex_latitude), np.max(self.aVertex_latitude))

                if not (-180 <= vertex_lon_range[0] <= 180 and -180 <= vertex_lon_range[1] <= 180):
                    logger.warning(f'Vertex longitudes outside valid range: {vertex_lon_range}')

                if not (-90 <= vertex_lat_range[0] <= 90 and -90 <= vertex_lat_range[1] <= 90):
                    logger.warning(f'Vertex latitudes outside valid range: {vertex_lat_range}')

            if not validation_passed:
                logger.error('Mesh topology rebuild failed validation')
                return None

            if iFlag_verbose:
                logger.info('Mesh topology successfully rebuilt')
                logger.info(f'Final mesh statistics:')
                logger.info(f'  - Unique vertices: {len(self.aVertex_longititude)}')
                logger.info(f'  - Mesh cells: {len(self.aCenter_longititude)}')
                logger.info(f'  - Max vertices per cell: {self.nVertex_max}')
                logger.info(f'  - Connectivity shape: {self.aConnectivity.shape}')
                logger.info(f'  - Vertex longitude range: [{np.min(self.aVertex_longititude):.3f}, {np.max(self.aVertex_longititude):.3f}]')
                logger.info(f'  - Vertex latitude range: [{np.min(self.aVertex_latitude):.3f}, {np.max(self.aVertex_latitude):.3f}]')

            return xv, yv, connectivity

        except Exception as e:
            logger.error(f'Unexpected error in rebuild_mesh_topology: {str(e)}')
            logger.error(f'Traceback: {traceback.format_exc()}')
            return None

    def report_inputs(self, iFlag_show_gpu_info=False):
        """
        Print comprehensive input information including raster and mesh details.

        Args:
            iFlag_show_gpu_info (bool): If True, also print GPU/GeoVista information
        """
        self.print_raster_info()
        self.print_mesh_info()

        if iFlag_show_gpu_info:
            try:
                import geovista.report as gvreport
                print("\n" + "="*60)
                print("GPU/GeoVista Information:")
                print("="*60)
                print(gvreport.Report())
            except ImportError:
                logger.warning("GeoVista not available for GPU info reporting")

    def report_outputs(self, sFilename_output=None):
        """
        Report output statistics.

        Args:
            sFilename_output (str, optional): Output file to report on
        """
        if sFilename_output and os.path.exists(sFilename_output):
            logger.info(f"Output file created: {sFilename_output}")
            logger.info(f"Output file size: {os.path.getsize(sFilename_output) / (1024*1024):.2f} MB")
        else:
            logger.warning("No output file information available")

    def print_raster_info(self):
        """
        Print detailed information about all input raster files.
        """
        print("\n" + "="*60)
        print(f"Input Raster Information ({len(self.aFilename_source_raster)} file(s)):")
        print("="*60)

        for idx, sFilename in enumerate(self.aFilename_source_raster, 1):
            print(f"\n[{idx}] {sFilename}")
            try:
                pRaster = sraster(sFilename)
                pRaster.read_metadata()
                pRaster.print_info()
            except Exception as e:
                logger.error(f"Error reading raster info: {e}")

    def print_mesh_info(self):
        """
        Print detailed mesh topology information.
        """
        if self.aCenter_longititude is None or len(self.aCenter_longititude) == 0:
            logger.warning("Mesh topology not yet built")
            return

        print("\n" + "="*60)
        print("Mesh Topology Information:")
        print("="*60)
        print(f"Number of mesh cells: {len(self.aCenter_longititude)}")
        print(f"Cell longitude range: {self.aCenter_longititude.min():.3f} to {self.aCenter_longititude.max():.3f}")
        print(f"Cell latitude range: {self.aCenter_latitude.min():.3f} to {self.aCenter_latitude.max():.3f}")
        print(f"Maximum vertices per cell: {self.nVertex_max}")

        if self.aVertex_longititude is not None:
            print(f"Total unique vertices: {len(self.aVertex_longititude)}")
        if self.aConnectivity is not None:
            print(f"Connectivity matrix shape: {self.aConnectivity.shape}")

        # Display area statistics if available
        if self.dArea_min is not None and self.dArea_max is not None:
            print(f"\nCell Area Statistics:")
            print(f"  Min area: {self.dArea_min:.6f}")
            print(f"  Max area: {self.dArea_max:.6f}")
            print(f"  Mean area: {self.dArea_mean:.6f}")
            print(f"  Area range ratio: {self.dArea_max/self.dArea_min:.2f}x" if self.dArea_min > 0 else "  Area range ratio: N/A")

        print("="*60)

    def run_remap(self, sFilename_target_mesh_out = None,
                    sFilename_source_mesh_in = None,
                    aFilename_source_raster_in = None,
                  iFlag_stat_in = 1,
                  iFlag_remap_method_in = 1,
                  iFlag_save_clipped_raster_in=0,
                  sFolder_raster_out_in = None,
                  iFlag_verbose=False):
        """
        Perform zonal statistics by clipping raster data to mesh polygons.

        This method delegates to the extract module for implementation.
        """
        if aFilename_source_raster_in is None:
            aFilename_source_raster = self.aFilename_source_raster
        else:
            aFilename_source_raster = aFilename_source_raster_in

        if sFilename_source_mesh_in is None:
            sFilename_source_mesh = self.sFilename_source_mesh
        else:
            sFilename_source_mesh = sFilename_source_mesh_in

        if sFilename_target_mesh_out is None:
            sFilename_target_mesh = self.sFilename_target_mesh
        else:
            sFilename_target_mesh = sFilename_target_mesh_out
            self.sFilename_target_mesh = sFilename_target_mesh_out
        return extract.run_remap(
             sFilename_target_mesh,
               sFilename_source_mesh,
               aFilename_source_raster,
             self.dArea_mean,
             iFlag_remap_method_in = iFlag_remap_method_in,
            iFlag_stat_in = iFlag_stat_in,
              iFlag_save_clipped_raster_in=iFlag_save_clipped_raster_in,
              sFolder_raster_out_in=sFolder_raster_out_in,
              iFlag_verbose=iFlag_verbose  )

    def visualize_source_mesh(self,
                              sFilename_out=None,
                              dLongitude_focus_in=0.0,
                              dLatitude_focus_in=0.0,
                              dZoom_factor=0.7,
                              iFlag_show_coastlines=True,
                              iFlag_show_graticule=True,
                              iFlag_verbose=False):
        """
        Visualize the source mesh topology using GeoVista 3D globe rendering.

        Creates an interactive or saved 3D visualization of the unstructured mesh
        with proper geographic context including coastlines and coordinate grid.

        Args:
            sFilename_out (str, optional): Output screenshot file path.
                If None, displays interactive viewer. Supports formats: .png, .jpg, .svg
            dLongitude_focus_in (float, optional): Camera focal point longitude in degrees.
                Valid range: -180 to 180. Default is 0.0 (prime meridian).
            dLatitude_focus_in (float, optional): Camera focal point latitude in degrees.
                Valid range: -90 to 90. Default is 0.0 (equator).
            dZoom_factor (float, optional): Camera zoom level.
                Higher values zoom in. Default is 0.7.
            iFlag_show_coastlines (bool, optional): Show coastline overlay.
                Default is True.
            iFlag_show_graticule (bool, optional): Show coordinate grid with labels.
                Default is True.
            iFlag_verbose (bool, optional): If True, print detailed progress messages.
                If False, only print error messages. Default is False.

        Returns:
            bool: True if visualization successful, False otherwise

        Note:
            - Requires 'geovista' package: pip install geovista
            - Interactive mode requires display environment
            - Mesh topology must be built before visualization (call rebuild_mesh_topology first)
        """
        return _visual.visualize_source_mesh(
            self, sFilename_out, dLongitude_focus_in, dLatitude_focus_in,
            dZoom_factor, iFlag_show_coastlines, iFlag_show_graticule, iFlag_verbose
        )

    def visualize_raster(self):
        """
        Visualize source raster data using GeoVista.

        Note:
            Not yet implemented. Placeholder for future raster visualization.
        """
        return _visual.visualize_raster(self)

    def visualize_target_mesh(self, sVariable_in=None,
                               sUnit_in=None,
                               sFilename_out=None,
                               dLongitude_focus_in=0.0,
                               dLatitude_focus_in=0.0,
                               dZoom_factor=0.7,
                               iFlag_show_coastlines=True,
                               iFlag_show_graticule=True,
                               sColormap='viridis',
                               iFlag_create_animation=False,
                               iAnimation_frames=360,
                               dAnimation_speed=1.0,
                               sAnimation_format='mp4',
                               iFlag_verbose=False):
        """
        Visualize the target mesh with computed zonal statistics using GeoVista 3D rendering.

        Creates an interactive or saved 3D visualization of the mesh with cells colored
        by computed statistics (mean, min, max, std) from raster processing. Can also
        create rotating animations by generating multiple frames.

        Args:
            sVariable_in (str): Variable field name to visualize.
                Common values: 'mean', 'min', 'max', 'std', 'area'
            sUnit_in (str, optional): Unit label for the colorbar (e.g., 'mm', 'kg/m²').
                Default is empty string.
            sFilename_out (str, optional): Output screenshot file path.
                If None, displays interactive viewer. Supports: .png, .jpg, .svg
                For animations, this becomes the base filename (e.g., 'animation.mp4')
            dLongitude_focus_in (float, optional): Camera focal point longitude in degrees.
                Valid range: -180 to 180. Default is 0.0. For animations, this is the starting longitude.
            dLatitude_focus_in (float, optional): Camera focal point latitude in degrees.
                Valid range: -90 to 90. Default is 0.0.
            dZoom_factor (float, optional): Camera zoom level.
                Higher values zoom in. Default is 0.75.
            iFlag_show_coastlines (bool, optional): Show coastline overlay.
                Default is True.
            iFlag_show_graticule (bool, optional): Show coordinate grid with labels.
                Default is True.
            sColormap (str, optional): Matplotlib colormap name.
                Default is 'viridis'. Examples: 'plasma', 'coolwarm', 'jet', 'RdYlBu'
            iFlag_create_animation (bool, optional): Create rotating animation.
                Default is False. When True, generates frames for 360° rotation.
            iAnimation_frames (int, optional): Number of frames for 360° rotation.
                Default is 36 (10° per frame). More frames = smoother animation.
            dAnimation_speed (float, optional): Animation speed in degrees per frame.
                Default is 10.0. Calculated as 360 / iAnimation_frames if not specified.
            sAnimation_format (str, optional): Animation output format.
                Default is 'mp4'. Supports: 'mp4', 'gif', 'avi'
            iFlag_verbose (bool, optional): If True, print detailed progress messages.
                If False, only print error messages. Default is False.

        Returns:
            bool: True if visualization successful, False otherwise

        Raises:
            ImportError: If geovista package is not installed
            ValueError: If target mesh file or required data is not available

        Note:
            - Requires 'geovista' package: pip install geovista
            - Target mesh file must exist (created by run_remap method)
            - Specified variable must exist as a field in the target mesh
            - Interactive mode requires display environment
            - Animation mode requires 'imageio' package for video creation: pip install imageio[ffmpeg]
        """
        return _visual.visualize_target_mesh(
            self, sVariable_in, sUnit_in, sFilename_out, dLongitude_focus_in, dLatitude_focus_in,
            dZoom_factor, iFlag_show_coastlines, iFlag_show_graticule, sColormap,
            iFlag_create_animation, iAnimation_frames, dAnimation_speed, sAnimation_format, iFlag_verbose
        )

    def _create_rotation_animation(self, plotter, sFilename_out, dLongitude_start, dLatitude_focus,
                                    iAnimation_frames, dAnimation_speed, sAnimation_format, iFlag_verbose=False):
        """
        Create a rotating animation of the 3D globe visualization.

        This method delegates to the _visual module for implementation.
        """
        return _visual._create_rotation_animation(
            self, plotter, sFilename_out, dLongitude_start, dLatitude_focus,
            iAnimation_frames, dAnimation_speed, sAnimation_format, iFlag_verbose
        )

    def cleanup(self):
        """
        Cleanup method to release spatial reference objects and other resources.
        """
        try:
            if hasattr(self, 'pSpatialRef') and self.pSpatialRef is not None:
                self.pSpatialRef = None
                logger.debug('Spatial reference object cleaned up successfully')
        except Exception as e:
            logger.warning(f'Error during cleanup of spatial reference: {e}')
