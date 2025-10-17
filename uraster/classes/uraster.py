# Define a class named 'uraster'
import os
import logging
import traceback
import signal
import sys
import time
import numpy as np
from osgeo import gdal, ogr, osr
from multiprocessing import Pool, cpu_count
from concurrent.futures import ThreadPoolExecutor, as_completed

from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_extension
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
from pyearth.gis.geometry.split_polygon_cross_idl import split_polygon_cross_idl
from pyearth.gis.geometry.extract_unique_vertices_and_connectivity import extract_unique_vertices_and_connectivity
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename, get_vector_format_from_filename

from .sraster import sraster

# Set up logging for crash detection
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)
crs = "EPSG:4326"
def signal_handler(signum, frame):
    """Handle system signals (like SIGSEGV) for crash detection"""
    logger.error(f"Received signal {signum}. Process may have crashed.")
    logger.error(f"Current frame: {frame}")
    traceback.print_stack(frame)
    sys.exit(1)

# Register signal handlers for common crash signals
signal.signal(signal.SIGSEGV, signal_handler)  # Segmentation fault
signal.signal(signal.SIGABRT, signal_handler)  # Abort
signal.signal(signal.SIGFPE, signal_handler)   # Floating point exception

# Constants for processing thresholds
IDL_LONGITUDE_THRESHOLD = 100  # Degrees - threshold for detecting IDL crossing
MEMORY_WARNING_THRESHOLD = 80  # Percent - warn when memory usage exceeds this
WARP_TIMEOUT_SECONDS = 30      # Seconds - timeout for GDAL Warp operations
PROGRESS_REPORT_INTERVAL = 100 # Report progress every N features

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

    def setup(self):
        """
        Initialize and validate the uraster configuration.
        Checks raster files and mesh file for existence and validity.

        Returns:
            bool: True if setup successful, False otherwise
        """
        raster_check = self.check_raster_files()
        mesh_check = self.check_mesh_file()

        return raster_check is not None and mesh_check is not None

    def check_raster_files(self, aFilename_source_raster_in=None):
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

        logger.info('All raster files exist and are readable')

        # Phase 2: Process and convert rasters to WGS84
        aFilename_source_raster_out = []

        # Create WGS84 spatial reference for comparison
        try:
            pSpatialRef_wgs84 = osr.SpatialReference()
            pSpatialRef_wgs84.ImportFromEPSG(4326)
            wkt_wgs84 = pSpatialRef_wgs84.ExportToWkt()
        except Exception as e:
            logger.error(f'Failed to create WGS84 spatial reference: {e}')
            return None

        # Process each raster file
        for idx, sFilename_raster_in in enumerate(aFilename_source_raster, 1):
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
                    logger.info(f'  ✓ Already in WGS84 (EPSG:4326)')
                    aFilename_source_raster_out.append(sFilename_raster_in)
                else:
                    # Convert to WGS84
                    logger.info(f'  → Converting to WGS84 from {pRaster.pSpatialRef.GetName() if pRaster.pSpatialRef else "unknown CRS"}')
                    try:
                        pRaster_wgs84 = pRaster.convert_to_wgs84()

                        if pRaster_wgs84 is None or not hasattr(pRaster_wgs84, 'sFilename'):
                            logger.error(f'Conversion to WGS84 failed: {sFilename_raster_in}')
                            return None

                        if not os.path.exists(pRaster_wgs84.sFilename):
                            logger.error(f'Converted WGS84 file not found: {pRaster_wgs84.sFilename}')
                            return None

                        logger.info(f'  ✓ Converted to: {pRaster_wgs84.sFilename}')
                        aFilename_source_raster_out.append(pRaster_wgs84.sFilename)

                    except Exception as e:
                        logger.error(f'Error during WGS84 conversion: {sFilename_raster_in}: {e}')
                        logger.error(f'Traceback: {traceback.format_exc()}')
                        return None

                # Log raster summary
                logger.debug(f'  - Dimensions: {pRaster.nrow} x {pRaster.ncolumn} pixels')
                logger.debug(f'  - Data type: {pRaster.eType}')
                if hasattr(pRaster, 'dNoData'):
                    logger.debug(f'  - NoData value: {pRaster.dNoData}')

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

        logger.info(f'Successfully validated and prepared {len(aFilename_source_raster_out)} raster file(s)')

        return aFilename_source_raster_out

    def check_mesh_file(self):
        """
        Check if the source mesh file exists and build its topology.

        Returns:
            tuple or None: (vertices_lon, vertices_lat, connectivity) if successful, None otherwise
        """
        if not self.sFilename_source_mesh:
            logger.error("No source mesh filename provided")
            return None

        if not os.path.exists(self.sFilename_source_mesh):
            logger.error(f"Source mesh file does not exist: {self.sFilename_source_mesh}")
            return None

        return self.rebuild_mesh_topology()

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

    def _determine_optimal_resampling(self, dPixelWidth, dPixelHeight):
        """
        Determine optimal resampling method based on mesh and raster resolution comparison.

        Compares the characteristic mesh cell size with raster resolution to decide
        whether to use nearest neighbor (when raster is much finer) or weighted
        averaging (when mesh and raster resolutions are comparable).

        Args:
            dPixelWidth (float): Raster pixel width in degrees
            dPixelHeight (float): Raster pixel height in degrees (absolute value)

        Returns:
            tuple: (resample_method_string, resample_method_code)
                - resample_method_string: 'nearest', 'bilinear', 'cubic', 'average', etc.
                - resample_method_code: integer code (1, 2, or 3)
        """
        if self.dArea_mean is None or self.dArea_mean <= 0:
            logger.warning("Mesh area statistics not available, defaulting to nearest neighbor")
            return 'near', 1

        # Estimate characteristic mesh cell dimension (approximate square root of mean area)
        dMesh_characteristic_size = np.sqrt(self.dArea_mean)

        # Use the coarser of the two raster dimensions
        dRaster_resolution = max(abs(dPixelWidth), abs(dPixelHeight))

        # Calculate resolution ratio (mesh size / raster size)
        dResolution_ratio = dMesh_characteristic_size / dRaster_resolution

        logger.info("="*60)
        logger.info("Resolution Comparison Analysis:")
        logger.info(f"  Raster resolution: {dRaster_resolution:.6f} degrees ({dRaster_resolution*111:.2f} km at equator)")
        logger.info(f"  Mean mesh cell size: {dMesh_characteristic_size:.6f} degrees ({dMesh_characteristic_size*111:.2f} km at equator)")
        logger.info(f"  Resolution ratio (mesh/raster): {dResolution_ratio:.2f}")
        logger.info(f"  Threshold for weighted averaging: {self.dResolution_ratio_threshold:.2f}")

        # Decision logic
        if dResolution_ratio < self.dResolution_ratio_threshold:
            # Mesh cells are comparable to or smaller than raster resolution
            # Use weighted averaging to properly capture sub-pixel variations
            recommended_method = 'average'
            recommended_code = 3
            logger.warning(f"Mesh resolution is close to raster resolution (ratio: {dResolution_ratio:.2f})")
            logger.warning(f"Switching to WEIGHTED AVERAGING (average) for accuracy")
            logger.warning("Consider using higher resolution raster data for better results")
        else:
            # Raster is much finer than mesh - nearest neighbor is appropriate
            recommended_method = 'near'
            recommended_code = 1
            logger.info(f"Raster is significantly finer than mesh (ratio: {dResolution_ratio:.2f})")
            logger.info(f"Using NEAREST NEIGHBOR resampling (sufficient for this resolution ratio)")

        logger.info("="*60)

        return recommended_method, recommended_code


    def rebuild_mesh_topology(self):
        """
        Rebuild mesh topology from source mesh file by extracting vertices,
        connectivity, and centroids for unstructured mesh processing.

        Returns:
            tuple: (vertices_longitude, vertices_latitude, connectivity) or None on failure
        """

        try:
            # Open the input data source
            pDataset = ogr.Open(self.sFilename_source_mesh, 0)  # Read-only
            if pDataset is None:
                logger.error(f'Failed to open file: {self.sFilename_source_mesh}')
                return None

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

            aCellID= list(range(nFeatures))

            logger.info(f'Processing {nFeatures} features with {iFieldCount} fields')

            # Get the first field name (assuming it contains the data variable)
            sVariable = pLayerDefn.GetFieldDefn(0).GetName() if iFieldCount > 0 else None

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
                            logger.warning(f'Feature {iFeature_index} has invalid geometry, skipping')
                            iFeature_index += 1
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

                                    aCellID[iFeature_index] = int(field_value) if field_value is not None else iFeature_index
                                except (ValueError, TypeError) as e:
                                    logger.warning(f'Could not convert field value for feature {iFeature_index}: {e}')
                                    data_list.append(iFeature_index)
                            else:
                                data_list.append(iFeature_index)  # Use feature index as default
                        else:
                            logger.warning(f'Failed to extract coordinates from feature {iFeature_index}')
                            invalid_geometry_count += 1

                    except Exception as e:
                        logger.warning(f'Error processing feature {iFeature_index}: {str(e)}')
                        invalid_geometry_count += 1

                elif sGeometry_type in ['MULTIPOLYGON', 'POINT', 'LINESTRING']:
                    logger.warning(f'Geometry type {sGeometry_type} not supported in feature {iFeature_index}, skipping')
                    invalid_geometry_count += 1
                else:
                    logger.warning(f'Unknown geometry type {sGeometry_type} in feature {iFeature_index}, skipping')
                    invalid_geometry_count += 1

                iFeature_index += 1

            # Report processing statistics
            valid_features = len(lons_list)
            logger.info(f'Feature processing summary:')
            logger.info(f'  - Total features: {iFeature_index}')
            logger.info(f'  - Valid polygons: {valid_features}')
            logger.info(f'  - Invalid/skipped: {invalid_geometry_count}')
            logger.info(f'  - Success rate: {(valid_features/iFeature_index*100):.1f}%' if iFeature_index > 0 else '  - Success rate: 0%')

            # Clean up dataset
            pDataset = None

            if not lons_list:
                logger.error('No valid polygon features found in mesh file')
                return None

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
                logger.info('Extracting unique vertices and connectivity...')
                xv, yv, connectivity, vertex_to_index = extract_unique_vertices_and_connectivity(
                    lons_list, lats_list
                )

                if xv is None or yv is None or connectivity is None:
                    logger.error('Failed to extract unique vertices and connectivity')
                    return None

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
                logger.warning(f"aCellID length ({len(aCellID)}) doesn't match mesh cells ({len(cell_lons_1d)}), truncating to match")
                aCellID = aCellID[:len(cell_lons_1d)]

            self.aCellID = np.array(aCellID)
            # Calculate and store area statistics
            if area_list:
                area_array = np.array(area_list)
                valid_areas = area_array[area_array > 0]  # Exclude zero areas from statistics
                if len(valid_areas) > 0:
                    self.dArea_min = float(np.min(valid_areas))
                    self.dArea_max = float(np.max(valid_areas))
                    self.dArea_mean = float(np.mean(valid_areas))
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

    def run_remap(self, sFilename_vector_out,
                    sFilename_target_mesh_in = None,
                                aFilename_source_raster_in = None,
                  iFlag_stat_in = 1,
                  iFlag_save_clipped_raster_in=0,
                  sFolder_raster_out_in = None,
                  sFormat_in='GTiff'):
        """
        Perform zonal statistics by clipping raster data to mesh polygons.

        Main processing method that extracts raster values for each mesh cell polygon
        and computes statistics (mean, min, max, std, sum, count).

        Args:
            sFilename_vector_out (str): Output vector file path with computed statistics
            sFilename_target_mesh_in (str, optional): Input mesh polygon file.
                Defaults to configured target mesh.
            aFilename_source_raster_in (list, optional): List of source raster files.
                Defaults to configured source rasters.
            iFlag_stat_in (int, optional): Flag to compute statistics (1=yes, 0=no).
                Default is 1.
            iFlag_save_clipped_raster_in (int, optional): Flag to save clipped rasters (1=yes, 0=no).
                Default is 0.
            sFolder_raster_out_in (str, optional): Output folder for clipped rasters.
                Required if iFlag_save_clipped_raster_in=1.
            sFormat_in (str, optional): GDAL raster format. Default is 'GTiff'.

        Returns:
            None

        Note:
            - Handles IDL-crossing polygons automatically
            - Generates failure report for problematic features
            - Supports multiple input rasters (uses first in list)
        """

        if aFilename_source_raster_in is None:
            aFilename_source_raster = self.aFilename_source_raster
        else:
            aFilename_source_raster = aFilename_source_raster_in

        if sFilename_target_mesh_in is None:
            sFilename_target_mesh = self.sFilename_target_mesh
        else:
            sFilename_target_mesh = sFilename_target_mesh_in

        #check input files
        for sFilename_raster in aFilename_source_raster:
            if os.path.exists(sFilename_raster):
                pass
            else:
                print('The raster file does not exist!')
                return

        if os.path.exists(sFilename_target_mesh ):
            pass
        else:
            print('The vector mesh file does not exist!')
            return

        # Determine output vector format from filename extension
        sVectorDriverName = get_vector_driver_from_extension(sFilename_vector_out)
        pDriver_vector = ogr.GetDriverByName(sVectorDriverName)

        #check the input raster data format and decide gdal driver
        if sFormat_in is not None:
            sDriverName = sFormat_in
        else:
            sDriverName = 'GTiff'

        if os.path.exists(sFilename_vector_out):
            #remove the file using the vector driver
            pDriver_vector.DeleteDataSource(sFilename_vector_out)

        pDriver = gdal.GetDriverByName(sDriverName)

        #get the raster file extension
        sFilename_raster = aFilename_source_raster[0]  #just use the first raster to get the extension
        sExtension = os.path.splitext(sFilename_raster)[1]
        sName = os.path.basename(sFilename_raster)
        sRasterName_no_extension = os.path.splitext(sName)[0]

        #use the sraster class the check the raster
        dX_left = -180.0
        dX_right = 180.0
        dY_top = 90.0
        dY_bot = -90.0
        dPixelWidth = None
        pPixelHeight = None
        #narrow the range to speed up the processing
        #also get the highest resolution
        for sFilename_raster in aFilename_source_raster:
            sRaster = sraster(sFilename_raster)
            sRaster.read_metadata()
            eType = sRaster.eType
            dX_left = min(dX_left, sRaster.dLongitude_left)
            dX_right = max(dX_right, sRaster.dLongitude_right)
            dY_top = max(dY_top, sRaster.dLatitude_top)
            dY_bot = min(dY_bot, sRaster.dLatitude_bottom)
            dMissing_value = sRaster.dNoData
            if dPixelWidth is None or sRaster.pTransform[1] < dPixelWidth:
                dPixelWidth = sRaster.pTransform[1]

            if pPixelHeight is None or sRaster.pTransform[5] > pPixelHeight: #is pPixelHeight is negative, then it should be less than
                pPixelHeight = sRaster.pTransform[5]

        # Determine optimal resampling method based on resolution comparison
        # This will override iFlag_remap_method if mesh resolution is too coarse
        sRemap_method_auto, iRemap_method_auto = self._determine_optimal_resampling(dPixelWidth, abs(pPixelHeight))

        # Use automatically determined method if it's more conservative than user setting
        # Priority: weighted averaging > nearest neighbor
        if iRemap_method_auto == 3 and self.iFlag_remap_method != 3:
            logger.warning(f"Overriding user remap method ({self.iFlag_remap_method}) with automatic selection (3 - weighted average)")
            logger.warning("This is necessary due to mesh/raster resolution compatibility")
            sRemap_method = sRemap_method_auto
            iFlag_remap_method_used = iRemap_method_auto
        else:
            # Use user's preferred method
            if self.iFlag_remap_method == 1:
                sRemap_method = 'near'
            elif self.iFlag_remap_method == 2:
                sRemap_method = 'near'
            elif self.iFlag_remap_method == 3:
                sRemap_method = 'average'
            iFlag_remap_method_used = self.iFlag_remap_method
            logger.info(f"Using user-specified remap method: {sRemap_method}")


        #get the spatial reference of the mesh vector file
        pDataset_mesh = ogr.Open( sFilename_target_mesh )
        pLayer_mesh = pDataset_mesh.GetLayer(0)
        nFeature = pLayer_mesh.GetFeatureCount()
        pSpatialRef_target = pLayer_mesh.GetSpatialRef()

        pSpatialRef_target_wkt = pSpatialRef_target.ExportToWkt()

        #check whether the polygon has only one or more features
        if nFeature > 1:
            pass
        else:
            print('The polygon file has only one polygons!')
            return

        options = ['COMPRESS=DEFLATE', 'PREDICTOR=2']
        sDriverName = 'MEM'

        #create a polygon feature to save the output
        pDataset_out = pDriver_vector.CreateDataSource(sFilename_vector_out)
        pLayer_out = pDataset_out.CreateLayer('cell', pSpatialRef_target, ogr.wkbPolygon)
        pLayer_defn_out = pLayer_out.GetLayerDefn()
        pFeature_out = ogr.Feature(pLayer_defn_out)

        #add id, area and mean, min, max, std of the raster
        pLayer_out.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
        #define a field
        pField = ogr.FieldDefn('area', ogr.OFTReal)
        pField.SetWidth(32)
        pField.SetPrecision(2)
        pLayer_out.CreateField(pField)

        #in the future, we will also copy other attributes from the input geojson file

        if iFlag_stat_in == 1:
            pLayer_out.CreateField(ogr.FieldDefn('mean', ogr.OFTReal))
            pLayer_out.CreateField(ogr.FieldDefn('min', ogr.OFTReal))
            pLayer_out.CreateField(ogr.FieldDefn('max', ogr.OFTReal))
            pLayer_out.CreateField(ogr.FieldDefn('std', ogr.OFTReal))
        else:
            #treat the raster as categorical?
            pass

        # Pre-compute GDAL options to avoid repeated object creation
        # sRemap_method was already determined above based on resolution comparison
        gdal_warp_options_base = {
            'cropToCutline': True,
            'xRes': dPixelWidth,
            'yRes': abs(pPixelHeight),
            'dstSRS': pSpatialRef_target,
            'format': 'MEM',
            'resampleAlg': sRemap_method,
            'srcSRS': 'EPSG:4326',  # Explicitly set source CRS
        }

        # Enable GDAL multi-threading
        gdal.SetConfigOption('GDAL_NUM_THREADS', str(cpu_count()))

        # Batch processing variables
        i = 1

        # Pre-fetch all features for potential parallel processing
        aFeatures = []
        pLayer_mesh.ResetReading()
        pFeature_mesh = pLayer_mesh.GetNextFeature()
        while pFeature_mesh is not None:
            # Check if feature crosses the international date line
            pPolygon = pFeature_mesh.GetGeometryRef()
            aCoord = get_geometry_coordinates(pPolygon)
            dLon_min = np.min(aCoord[:,0])
            dLon_max = np.max(aCoord[:,0])

            if dLon_max - dLon_min > IDL_LONGITUDE_THRESHOLD:  # Use constant
                if not pPolygon.IsValid():
                    logger.warning(f'Feature {i} has invalid geometry crossing IDL')

                # Create multipolygon to handle IDL crossing
                logger.info(f'Feature {i} crosses the international date line, splitting into multiple parts.')
                pMultipolygon = ogr.Geometry(ogr.wkbMultiPolygon)
                aCoord_gcs_split = split_polygon_cross_idl(aCoord)

                for aCoord_gcs in aCoord_gcs_split:
                    #create a polygon (not just ring) and add it to the multipolygon
                    ring = ogr.Geometry(ogr.wkbLinearRing)
                    for iCoord in range(aCoord_gcs.shape[0]):
                        ring.AddPoint(aCoord_gcs[iCoord, 0], aCoord_gcs[iCoord, 1])

                    ring.CloseRings()
                    # Create polygon from ring
                    polygon_part = ogr.Geometry(ogr.wkbPolygon)
                    polygon_part.AddGeometry(ring)
                    pMultipolygon.AddGeometry(polygon_part)  # Add polygon, not ring

                #create a feature from the multipolygon and add it to the list
                pFeature_new = ogr.Feature(pLayer_mesh.GetLayerDefn())
                pFeature_new.SetGeometry(pMultipolygon)
                aFeatures.append(pFeature_new)
                pFeature_mesh = pLayer_mesh.GetNextFeature()
            else:
                aFeatures.append(pFeature_mesh.Clone())
                pFeature_mesh = pLayer_mesh.GetNextFeature()

            i += 1

        print(f"Processing {len(aFeatures)} features...")
        start_time = time.time()

        # Add crash tracking variables
        failed_features = []
        successful_features = 0

        #reset i
        i = 1

        # Process features
        for idx, pFeature_mesh in enumerate(aFeatures):
            sClip = f"{i:08d}"  # f-string is faster than format
            pPolygon = pFeature_mesh.GetGeometryRef()
            # Fast envelope check first
            minX, maxX, minY, maxY = pPolygon.GetEnvelope()
            if (minX > dX_right or maxX < dX_left or
                minY > dY_top or maxY < dY_bot or
                not pPolygon or pPolygon.IsEmpty() or not pPolygon.IsValid()):
                i += 1
                continue

            # Process valid polygons with comprehensive error handling
            feature_start_time = time.time()
            try:
                # Monitor memory usage (optional - requires psutil)
                try:
                    import psutil
                    process = psutil.Process()
                    memory_percent = process.memory_percent()
                    if memory_percent > MEMORY_WARNING_THRESHOLD:  # Use constant
                        logger.warning(f"High memory usage: {memory_percent:.1f}% at feature {i}")
                except ImportError:
                    pass  # psutil not available

                # Handle polygon vs multipolygon differently
                geometry_type = pPolygon.GetGeometryType()
                sGeometry_type = pPolygon.GetGeometryName()
                geometry_type_name = self._get_geometry_type_name(geometry_type)
                logger.debug(f"Feature {i} geometry type: {geometry_type} ({geometry_type_name})")

                if sGeometry_type == "POLYGON":
                    # Simple polygon - process normally
                    pDataset_clip_warped, aData_clip, newGeoTransform = self._process_single_polygon(
                        pPolygon, aFilename_source_raster, gdal_warp_options_base, i)

                    if pDataset_clip_warped is None:
                        error_msg = f"Failed to process single polygon for feature {i}"
                        logger.error(error_msg)
                        failed_features.append({'feature_id': i, 'error': error_msg, 'envelope': (minX, maxX, minY, maxY)})
                        i += 1
                        continue

                elif sGeometry_type == "MULTIPOLYGON":
                    # Multipolygon (IDL crossing) - process each part separately and merge
                    logger.info(f"Processing IDL-crossing multipolygon for feature {i}")
                    merged_data, merged_transform = self._process_multipolygon_idl(
                        pPolygon, aFilename_source_raster, gdal_warp_options_base, i, dMissing_value)

                    if merged_data is None:
                        error_msg = f"Failed to process IDL-crossing multipolygon for feature {i}"
                        logger.error(error_msg)
                        failed_features.append({'feature_id': i, 'error': error_msg, 'envelope': (minX, maxX, minY, maxY)})
                        i += 1
                        continue

                    aData_clip = merged_data
                    newGeoTransform = merged_transform
                    pDataset_clip_warped = True  # Flag to indicate successful processing

                else:
                    error_msg = f"Unsupported geometry type for feature {i}: {geometry_type} ({geometry_type_name})"
                    logger.error(error_msg)
                    failed_features.append({'feature_id': i, 'error': error_msg, 'envelope': (minX, maxX, minY, maxY)})
                    i += 1
                    continue

                logger.debug(f"Successfully processed GDAL Warp for feature {i} in {time.time() - feature_start_time:.2f}s")
                successful_features += 1

            except TimeoutError as e:
                logger.error(str(e))
                failed_features.append({'feature_id': i, 'error': str(e), 'envelope': (minX, maxX, minY, maxY)})
                i += 1
                continue

            except MemoryError as e:
                error_msg = f"Memory error during GDAL Warp for feature {i}: {str(e)}"
                logger.error(error_msg)
                failed_features.append({'feature_id': i, 'error': error_msg, 'envelope': (minX, maxX, minY, maxY)})
                # Force garbage collection
                import gc
                gc.collect()
                i += 1
                continue

            except Exception as e:
                error_msg = f"Unexpected exception during GDAL Warp for feature {i}: {str(e)}"
                logger.error(error_msg)
                logger.error(f"Polygon envelope: {minX}, {maxX}, {minY}, {maxY}")
                logger.error(f"Exception type: {type(e).__name__}")
                logger.error(f"Traceback: {traceback.format_exc()}")
                # Log polygon geometry for debugging (truncated)
                try:
                    pPolygonWKT = pPolygon.ExportToWkt()
                    logger.debug(f"Polygon WKT (first 200 chars): {pPolygonWKT[:200]}...")
                except:
                    logger.debug("Could not log polygon WKT")
                failed_features.append({'feature_id': i, 'error': error_msg, 'envelope': (minX, maxX, minY, maxY)})
                i += 1
                continue

            # Optimize data type conversion and missing value handling
            if aData_clip is not None:
                # Use numpy's faster operations
                aData_clip = aData_clip.astype(np.int32, copy=False)  # Avoid unnecessary copy
                np.place(aData_clip, aData_clip == dMissing_value, -9999)  # Faster than boolean indexing

                if iFlag_save_clipped_raster_in == 1 :
                    if geometry_type == ogr.wkbPolygon : #we cannot save the multipolygon directly
                        sFilename_raster_out = os.path.join(sFolder_raster_out_in, f'{sRasterName_no_extension}_clip_{sClip}{sExtension}')
                        # Delete existing file if it exists (GDAL Create() doesn't overwrite)
                        if os.path.exists(sFilename_raster_out):
                            try:
                                pDriver.Delete(sFilename_raster_out)
                            except:
                                os.remove(sFilename_raster_out)  # Fallback if GDAL delete fails
                        iNewWidth = aData_clip.shape[1]
                        iNewHeigh = aData_clip.shape[0]
                        pDataset_clip = pDriver.Create(sFilename_raster_out, iNewWidth, iNewHeigh, 1, eType , options= options)
                        pDataset_clip.SetGeoTransform( newGeoTransform )
                        pDataset_clip.SetProjection( pSpatialRef_target_wkt)
                        #set the no data value
                        pBand = pDataset_clip.GetRasterBand(1)
                        pBand.SetNoDataValue(-9999)
                        pBand.WriteArray(aData_clip)
                        pBand.FlushCache()  # Corrected method name to FlushCache()
                        pDataset_clip.FlushCache()
                        pDataset_clip = None

                aCoords_gcs = get_geometry_coordinates(pPolygon)
                dArea = calculate_polygon_area(aCoords_gcs[:,0], aCoords_gcs[:,1])
                #create a polygon feature to save the output
                pFeature_out.SetGeometry(pPolygon.Clone())
                pFeature_out.SetField('id', i)
                pFeature_out.SetField('area', dArea)

                if iFlag_stat_in == 1:
                    # Check if aData_clip is already 1D valid data (from IDL merge) or 2D raster
                    if aData_clip.ndim == 1:
                        # Already 1D array of valid data from IDL merge
                        valid_data = aData_clip[aData_clip != -9999]  # Remove any remaining -9999 values
                    else:
                        # Standard 2D raster - extract valid data
                        valid_mask = aData_clip != -9999
                        valid_data = aData_clip[valid_mask]

                    if valid_data.size > 0:
                        # Calculate all statistics in one pass for efficiency
                        valid_data_float = valid_data.astype(np.float64, copy=False)
                        pFeature_out.SetField('mean', float(np.mean(valid_data_float)))
                        pFeature_out.SetField('min', float(np.min(valid_data_float)))
                        pFeature_out.SetField('max', float(np.max(valid_data_float)))
                        pFeature_out.SetField('std', float(np.std(valid_data_float)))

                pLayer_out.CreateFeature(pFeature_out)

                # Explicit cleanup - handle both single polygon dataset and multipolygon flag
                if isinstance(pDataset_clip_warped, bool):
                    # Multipolygon case - cleanup was already handled in helper method
                    pass
                else:
                    # Single polygon case - cleanup the dataset
                    pDataset_clip_warped = None

                # Progress reporting
                if i % PROGRESS_REPORT_INTERVAL == 0:
                    elapsed = time.time() - start_time
                    rate = i / elapsed
                    eta = (len(aFeatures) - i) / rate if rate > 0 else 0
                    print(f"Processed {i}/{len(aFeatures)} features ({rate:.2f} features/sec, ETA: {eta:.0f}s)")

            i += 1


        pDataset_out.FlushCache()  # Flush once after all features are added
        pDataset_out = None        # Close the dataset
        pDataset_mesh = None

        # Report processing summary
        total_time = time.time() - start_time
        logger.info(f"Processing completed in {total_time:.2f} seconds")
        logger.info(f"Successfully processed: {successful_features} features")
        logger.info(f"Failed features: {len(failed_features)}")

        if failed_features:
            logger.warning("Failed features summary:")
            for failed in failed_features[:10]:  # Show first 10 failures
                logger.warning(f"  Feature {failed['feature_id']}: {failed['error']}")
            if len(failed_features) > 10:
                logger.warning(f"  ... and {len(failed_features) - 10} more failures")

        # Save failure report to file
        if failed_features:
            # Generate failure report filename by replacing extension with '_failures.log'
            base_name = os.path.splitext(sFilename_vector_out)[0]
            failure_report_file = f"{base_name}_failures.log"
            try:
                with open(failure_report_file, 'w') as f:
                    f.write(f"Processing failure report - {time.ctime()}\n")
                    f.write(f"Total features processed: {len(aFeatures)}\n")
                    f.write(f"Successful: {successful_features}\n")
                    f.write(f"Failed: {len(failed_features)}\n\n")
                    for failed in failed_features:
                        f.write(f"Feature {failed['feature_id']}: {failed['error']}\n")
                        f.write(f"  Envelope: {failed['envelope']}\n\n")
                logger.info(f"Failure report saved to: {failure_report_file}")
            except Exception as e:
                logger.error(f"Could not save failure report: {e}")

        return

    def visualize_source_mesh(self,
                              sFilename_out=None,
                              dLongitude_focus_in=0.0,
                              dLatitude_focus_in=0.0,
                              dZoom_factor=1.4,
                              iFlag_show_coastlines=True,
                              iFlag_show_graticule=True):
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
                Higher values zoom in. Default is 1.4.
            iFlag_show_coastlines (bool, optional): Show coastline overlay.
                Default is True.
            iFlag_show_graticule (bool, optional): Show coordinate grid with labels.
                Default is True.

        Returns:
            bool: True if visualization successful, False otherwise

        Note:
            - Requires 'geovista' package: pip install geovista
            - Interactive mode requires display environment
            - Mesh topology must be built before visualization (call rebuild_mesh_topology first)
        """
        # Validate mesh data availability
        if self.aVertex_longititude is None or self.aVertex_latitude is None:
            logger.error('Mesh vertices not available. Build mesh topology first.')
            return False

        if self.aConnectivity is None:
            logger.error('Mesh connectivity not available. Build mesh topology first.')
            return False

        mesh.cell_data[name] = self.aCellID

        # Plot the mesh.
        plotter = gv.GeoPlotter()
        sargs = {"title": f"{name} / {sUnit}",
               "shadow": True,
                "title_font_size": 10,
                    "label_font_size": 10,
                        "fmt": "%.1f",
        }
        plotter.add_mesh(mesh, scalars=name, scalar_bar_args=sargs)
        focal_point = gv.geodesic.to_cartesian([dLongitude_focus], [dLatitude_focus])[0]
        camera_position = gv.geodesic.to_cartesian([dLongitude_focus], [dLatitude_focus], radius=gv.common.RADIUS *6)[0]
        plotter.camera.focal_point = focal_point
        plotter.camera.position = camera_position
        plotter.camera.zoom(1.4)
        plotter.add_coastlines()
        plotter.add_axes()
        plotter.add_graticule(show_labels=True)
        if len(self.aVertex_longititude) == 0 or len(self.aVertex_latitude) == 0:
            logger.error('Mesh vertices are empty.')
            return False

        if self.aConnectivity.size == 0:
            logger.error('Mesh connectivity is empty.')
            return False

        # Validate focus coordinates
        dLongitude_focus = dLongitude_focus_in if dLongitude_focus_in is not None else 0.0
        dLatitude_focus = dLatitude_focus_in if dLatitude_focus_in is not None else 0.0

        if not (-180 <= dLongitude_focus <= 180):
            logger.warning(f'Longitude focus {dLongitude_focus} out of range [-180, 180], clamping')
            dLongitude_focus = np.clip(dLongitude_focus, -180, 180)

        if not (-90 <= dLatitude_focus <= 90):
            logger.warning(f'Latitude focus {dLatitude_focus} out of range [-90, 90], clamping')
            dLatitude_focus = np.clip(dLatitude_focus, -90, 90)

        # Validate zoom factor
        if dZoom_factor <= 0:
            logger.warning(f'Invalid zoom factor {dZoom_factor}, using default 1.4')
            dZoom_factor = 1.4

        # Validate output file path if provided
        if sFilename_out is not None:
            if not isinstance(sFilename_out, str) or not sFilename_out.strip():
                logger.error('Output filename must be a non-empty string')
                return False

            # Check output directory exists
            output_dir = os.path.dirname(sFilename_out)
            if output_dir and not os.path.exists(output_dir):
                try:
                    os.makedirs(output_dir, exist_ok=True)
                    logger.info(f'Created output directory: {output_dir}')
                except Exception as e:
                    logger.error(f'Cannot create output directory {output_dir}: {e}')
                    return False

            # Check supported file extensions
            valid_extensions = ['.png', '.jpg', '.jpeg', '.svg', '.tif', '.tiff']
            file_ext = os.path.splitext(sFilename_out)[1].lower()
            if file_ext not in valid_extensions:
                logger.warning(f'File extension {file_ext} may not be supported. Recommended: .png, .jpg, .svg')

        # Import geovista with error handling
        try:
            import geovista as gv
            logger.info('GeoVista library imported successfully')
        except ImportError as e:
            logger.error('GeoVista library not available. Install with: pip install geovista')
            logger.error(f'Import error: {e}')
            return False

        try:
            logger.info('Creating mesh visualization...')
            logger.info(f'  - Vertices: {len(self.aVertex_longititude)}')
            logger.info(f'  - Connectivity shape: {self.aConnectivity.shape}')
            logger.info(f'  - Focus: ({dLongitude_focus:.2f}°, {dLatitude_focus:.2f}°)')
            logger.info(f'  - Zoom factor: {dZoom_factor}')

            # Prepare mesh metadata
            name = 'Mesh Cell ID'
            sUnit = ""

            # Create masked connectivity array (mask invalid indices)
            connectivity_masked = np.ma.masked_where(
                self.aConnectivity == -1,
                self.aConnectivity
            )

            # Validate connectivity indices
            valid_connectivity = self.aConnectivity[self.aConnectivity >= 0]
            if len(valid_connectivity) > 0:
                max_vertex_idx = len(self.aVertex_longititude) - 1
                if np.max(valid_connectivity) > max_vertex_idx:
                    logger.error(f'Connectivity contains invalid vertex index: max={np.max(valid_connectivity)}, vertices={len(self.aVertex_longititude)}')
                    return False

            # Transform to GeoVista unstructured mesh
            mesh = gv.Transform.from_unstructured(
                self.aVertex_longititude,
                self.aVertex_latitude,
                connectivity=connectivity_masked,
                crs=crs
            )

            logger.info(f'Created GeoVista mesh with {mesh.n_cells} cells and {mesh.n_points} points')

            # Create 3D plotter
            plotter = gv.GeoPlotter()

            # Configure scalar bar (colorbar) appearance
            sargs = {
                "title": f"{name} / {sUnit}" if sUnit else name,
                "shadow": True,
                "title_font_size": 10,
                "label_font_size": 10,
                "fmt": "%.0f",  # Integer formatting for cell IDs
                "n_labels": 5,
            }

            # Add mesh to plotter
            plotter.add_mesh(mesh, scalars=name, scalar_bar_args=sargs)

            # Configure camera position and focus
            try:
                focal_point = gv.geodesic.to_cartesian(
                    [dLongitude_focus],
                    [dLatitude_focus]
                )[0]

                camera_position = gv.geodesic.to_cartesian(
                    [dLongitude_focus],
                    [dLatitude_focus],
                    radius=gv.common.RADIUS * 6
                )[0]

                plotter.camera.focal_point = focal_point
                plotter.camera.position = camera_position
                plotter.camera.zoom(dZoom_factor)

                logger.debug(f'Camera configured: focal={focal_point}, position={camera_position}')
            except Exception as e:
                logger.warning(f'Error setting camera position: {e}. Using default view.')

            # Add geographic context
            if iFlag_show_coastlines:
                try:
                    plotter.add_coastlines()
                    logger.debug('Added coastlines overlay')
                except Exception as e:
                    logger.warning(f'Could not add coastlines: {e}')

            # Add coordinate axes
            try:
                plotter.add_axes()
                logger.debug('Added coordinate axes')
            except Exception as e:
                logger.warning(f'Could not add axes: {e}')

            # Add graticule (coordinate grid)
            if iFlag_show_graticule:
                try:
                    plotter.add_graticule(show_labels=True)
                    logger.debug('Added coordinate graticule with labels')
                except Exception as e:
                    logger.warning(f'Could not add graticule: {e}')

            # Output or display
            if sFilename_out is not None:
                # Save screenshot
                try:
                    plotter.screenshot(sFilename_out)
                    logger.info(f'✓ Visualization saved to: {sFilename_out}')

                    # Verify file was created
                    if os.path.exists(sFilename_out):
                        file_size = os.path.getsize(sFilename_out)
                        logger.info(f'  File size: {file_size / 1024:.1f} KB')
                    else:
                        logger.warning(f'Screenshot command executed but file not found: {sFilename_out}')

                    plotter.close()
                    return True

                except Exception as e:
                    logger.error(f'Failed to save screenshot: {e}')
                    logger.error(f'Traceback: {traceback.format_exc()}')
                    plotter.close()
                    return False
            else:
                # Interactive display
                try:
                    logger.info('Opening interactive visualization window...')
                    plotter.show()
                    return True
                except Exception as e:
                    logger.error(f'Failed to display interactive visualization: {e}')
                    logger.error(f'Ensure display environment is available (X11, Wayland, etc.)')
                    logger.error(f'Traceback: {traceback.format_exc()}')
                    plotter.close()
                    return False

        except ImportError as e:
            logger.error(f'Missing required GeoVista dependencies: {e}')
            return False

        except Exception as e:
            logger.error(f'Unexpected error during mesh visualization: {e}')
            logger.error(f'Error type: {type(e).__name__}')
            logger.error(f'Traceback: {traceback.format_exc()}')
            return False

    def visualize_raster(self):
        """
        Visualize source raster data using GeoVista.

        Note:
            Not yet implemented. Placeholder for future raster visualization.
        """
        #it is possible to visualize the raster using geovista as well
        #if more than one rasters is used, then we may overlap them one by one
        for idx, sFilename in enumerate(self.aFilename_source_raster, 1):
            print(f"\n[{idx}] {sFilename}")
            try:
                pRaster = sraster(sFilename)
                pRaster.create_raster_mesh()
                sFilename_raster_mesh = pRaster.sFilename_mesh

                #use this mesh to visualize, with the raster value as the cell data


            except Exception as e:
                logger.error(f"Error reading raster info: {e}")
        return

    def visualize_target_mesh(self, sVariable_in,
                               sUnit_in=None,
                               sFilename_out=None,
                               dLongitude_focus_in=0.0,
                               dLatitude_focus_in=0.0,
                               dZoom_factor=1.4,
                               iFlag_show_coastlines=True,
                               iFlag_show_graticule=True,
                               sColormap='viridis'):
        """
        Visualize the target mesh with computed zonal statistics using GeoVista 3D rendering.

        Creates an interactive or saved 3D visualization of the mesh with cells colored
        by computed statistics (mean, min, max, std) from raster processing.

        Args:
            sVariable_in (str): Variable field name to visualize.
                Common values: 'mean', 'min', 'max', 'std', 'area'
            sUnit_in (str, optional): Unit label for the colorbar (e.g., 'mm', 'kg/m²').
                Default is empty string.
            sFilename_out (str, optional): Output screenshot file path.
                If None, displays interactive viewer. Supports: .png, .jpg, .svg
            dLongitude_focus_in (float, optional): Camera focal point longitude in degrees.
                Valid range: -180 to 180. Default is 0.0.
            dLatitude_focus_in (float, optional): Camera focal point latitude in degrees.
                Valid range: -90 to 90. Default is 0.0.
            dZoom_factor (float, optional): Camera zoom level.
                Higher values zoom in. Default is 1.4.
            iFlag_show_coastlines (bool, optional): Show coastline overlay.
                Default is True.
            iFlag_show_graticule (bool, optional): Show coordinate grid with labels.
                Default is True.
            sColormap (str, optional): Matplotlib colormap name.
                Default is 'viridis'. Examples: 'plasma', 'coolwarm', 'jet', 'RdYlBu'

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
        """
        # Validate target mesh file
        if not self.sFilename_target_mesh:
            logger.error('No target mesh filename configured')
            return False

        if not os.path.exists(self.sFilename_target_mesh):
            logger.error(f'Target mesh file does not exist: {self.sFilename_target_mesh}')
            return False

        # Validate mesh topology data
        if self.aVertex_longititude is None or self.aVertex_latitude is None:
            logger.error('Mesh vertices not available. Build mesh topology first.')
            return False

        if self.aConnectivity is None:
            logger.error('Mesh connectivity not available. Build mesh topology first.')
            return False

        if len(self.aVertex_longititude) == 0 or len(self.aVertex_latitude) == 0:
            logger.error('Mesh vertices are empty.')
            return False

        # Validate variable name
        if not sVariable_in:
            logger.warning('No variable specified, defaulting to "mean"')
            sVariable = 'mean'
        else:
            sVariable = sVariable_in

        if not isinstance(sVariable, str):
            logger.error(f'Variable name must be a string, got {type(sVariable).__name__}')
            return False

        # Validate focus coordinates
        dLongitude_focus = dLongitude_focus_in if dLongitude_focus_in is not None else 0.0
        dLatitude_focus = dLatitude_focus_in if dLatitude_focus_in is not None else 0.0

        if not (-180 <= dLongitude_focus <= 180):
            logger.warning(f'Longitude focus {dLongitude_focus} out of range [-180, 180], clamping')
            dLongitude_focus = np.clip(dLongitude_focus, -180, 180)

        if not (-90 <= dLatitude_focus <= 90):
            logger.warning(f'Latitude focus {dLatitude_focus} out of range [-90, 90], clamping')
            dLatitude_focus = np.clip(dLatitude_focus, -90, 90)

        # Validate zoom factor
        if dZoom_factor <= 0:
            logger.warning(f'Invalid zoom factor {dZoom_factor}, using default 1.4')
            dZoom_factor = 1.4

        # Validate output file path if provided
        if sFilename_out is not None:
            if not isinstance(sFilename_out, str) or not sFilename_out.strip():
                logger.error('Output filename must be a non-empty string')
                return False

            # Check output directory exists
            output_dir = os.path.dirname(sFilename_out)
            if output_dir and not os.path.exists(output_dir):
                try:
                    os.makedirs(output_dir, exist_ok=True)
                    logger.info(f'Created output directory: {output_dir}')
                except Exception as e:
                    logger.error(f'Cannot create output directory {output_dir}: {e}')
                    return False

            # Check supported file extensions
            valid_extensions = ['.png', '.jpg', '.jpeg', '.svg', '.tif', '.tiff']
            file_ext = os.path.splitext(sFilename_out)[1].lower()
            if file_ext not in valid_extensions:
                logger.warning(f'File extension {file_ext} may not be supported. Recommended: .png, .jpg, .svg')

        # Import geovista with error handling
        try:
            import geovista as gv
            logger.info('GeoVista library imported successfully')
        except ImportError as e:
            logger.error('GeoVista library not available. Install with: pip install geovista')
            logger.error(f'Import error: {e}')
            return False

        try:
            logger.info(f'Loading target mesh data from: {self.sFilename_target_mesh}')

            # Open target mesh file
            pDataset = ogr.Open(self.sFilename_target_mesh, 0)  # Read-only
            if pDataset is None:
                logger.error(f'Failed to open target mesh file: {self.sFilename_target_mesh}')
                return False

            # Get first layer
            pLayer = pDataset.GetLayer(0)
            if pLayer is None:
                logger.error('Failed to get layer from target mesh dataset')
                pDataset = None
                return False

            # Get layer definition
            pLayerDefn = pLayer.GetLayerDefn()
            if pLayerDefn is None:
                logger.error('Failed to get layer definition from target mesh')
                pDataset = None
                return False

            # Get field information
            iFieldCount = pLayerDefn.GetFieldCount()
            nFeatures = pLayer.GetFeatureCount()
            self.nCell_target = nFeatures

            logger.info(f'Target mesh contains {nFeatures} features with {iFieldCount} fields')

            # Check if variable field exists
            field_names = [pLayerDefn.GetFieldDefn(i).GetName() for i in range(iFieldCount)]
            if sVariable not in field_names:
                logger.error(f'Variable "{sVariable}" not found in target mesh')
                logger.error(f'Available fields: {", ".join(field_names)}')
                pDataset = None
                return False

            logger.info(f'Extracting variable: {sVariable}')

            # Validate feature count matches source
            if self.nCell_source > 0 and self.nCell_target != self.nCell_source:
                logger.warning(f'Feature count mismatch: target={self.nCell_target}, source={self.nCell_source}')

            # Extract variable data from features
            data_list = []
            pLayer.ResetReading()
            pFeature = pLayer.GetNextFeature()
            feature_count = 0

            while pFeature is not None:
                pGeometry = pFeature.GetGeometryRef()
                if pGeometry is not None:
                    sGeometry_type = pGeometry.GetGeometryName()
                    if sGeometry_type == 'POLYGON':
                        try:
                            field_value = pFeature.GetField(sVariable)
                            if field_value is not None:
                                data_list.append(field_value)
                            else:
                                # Handle None values as NaN
                                data_list.append(np.nan)
                                logger.warning(f'Feature {feature_count} has None value for {sVariable}')
                        except Exception as e:
                            logger.warning(f'Error reading field {sVariable} from feature {feature_count}: {e}')
                            data_list.append(np.nan)
                    else:
                        logger.warning(f'Feature {feature_count} has unsupported geometry type: {sGeometry_type}')
                        data_list.append(np.nan)
                else:
                    logger.warning(f'Feature {feature_count} has no geometry')
                    data_list.append(np.nan)

                feature_count += 1
                pFeature = pLayer.GetNextFeature()

            # Close dataset
            pDataset = None

            if not data_list:
                logger.error('No data extracted from target mesh')
                return False

            # Convert to numpy array
            data = np.array(data_list, dtype=np.float64)

            # Validate data
            valid_data_mask = np.isfinite(data)
            valid_data_count = np.sum(valid_data_mask)

            if valid_data_count == 0:
                logger.error(f'All values for variable "{sVariable}" are invalid (NaN/Inf)')
                return False

            if valid_data_count < len(data):
                logger.warning(f'{len(data) - valid_data_count} of {len(data)} values are invalid')

            # Log data statistics
            valid_values = data[valid_data_mask]
            logger.info(f'Data statistics for "{sVariable}":')
            logger.info(f'  - Valid values: {valid_data_count}/{len(data)}')
            logger.info(f'  - Min: {np.min(valid_values):.4f}')
            logger.info(f'  - Max: {np.max(valid_values):.4f}')
            logger.info(f'  - Mean: {np.mean(valid_values):.4f}')
            logger.info(f'  - Std: {np.std(valid_values):.4f}')

            # Create masked connectivity array
            connectivity_masked = np.ma.masked_where(
                self.aConnectivity == -1,
                self.aConnectivity
            )

            # Validate connectivity indices
            valid_connectivity = self.aConnectivity[self.aConnectivity >= 0]
            if len(valid_connectivity) > 0:
                max_vertex_idx = len(self.aVertex_longititude) - 1
                if np.max(valid_connectivity) > max_vertex_idx:
                    logger.error(f'Connectivity contains invalid vertex index: max={np.max(valid_connectivity)}, vertices={len(self.aVertex_longititude)}')
                    return False

            # Transform to GeoVista unstructured mesh
            logger.info('Creating GeoVista mesh...')
            mesh = gv.Transform.from_unstructured(
                self.aVertex_longititude,
                self.aVertex_latitude,
                connectivity=connectivity_masked,
                crs=crs
            )

            logger.info(f'Created mesh with {mesh.n_cells} cells and {mesh.n_points} points')

            # Attach data to mesh
            name = sVariable.capitalize()
            sUnit = sUnit_in if sUnit_in is not None else ""
            mesh.cell_data[name] = data

            # Create 3D plotter
            plotter = gv.GeoPlotter()

            # Configure scalar bar (colorbar) appearance
            sargs = {
                "title": f"{name} / {sUnit}" if sUnit else name,
                "shadow": True,
                "title_font_size": 10,
                "label_font_size": 10,
                "fmt": "%.2f",  # Floating point formatting for statistics
                "n_labels": 5,
            }

            # Add mesh to plotter with colormap
            plotter.add_mesh(mesh, scalars=name, scalar_bar_args=sargs, cmap=sColormap)

            # Configure camera position and focus
            try:
                focal_point = gv.geodesic.to_cartesian(
                    [dLongitude_focus],
                    [dLatitude_focus]
                )[0]

                camera_position = gv.geodesic.to_cartesian(
                    [dLongitude_focus],
                    [dLatitude_focus],
                    radius=gv.common.RADIUS * 6
                )[0]

                plotter.camera.focal_point = focal_point
                plotter.camera.position = camera_position
                plotter.camera.zoom(dZoom_factor)

                logger.debug(f'Camera configured: focal={focal_point}, position={camera_position}')
            except Exception as e:
                logger.warning(f'Error setting camera position: {e}. Using default view.')

            # Add geographic context
            if iFlag_show_coastlines:
                try:
                    plotter.add_coastlines()
                    logger.debug('Added coastlines overlay')
                except Exception as e:
                    logger.warning(f'Could not add coastlines: {e}')

            # Add coordinate axes
            try:
                plotter.add_axes()
                logger.debug('Added coordinate axes')
            except Exception as e:
                logger.warning(f'Could not add axes: {e}')

            # Add graticule (coordinate grid)
            if iFlag_show_graticule:
                try:
                    plotter.add_graticule(show_labels=True)
                    logger.debug('Added coordinate graticule with labels')
                except Exception as e:
                    logger.warning(f'Could not add graticule: {e}')

            # Output or display
            if sFilename_out is not None:
                # Save screenshot
                try:
                    plotter.screenshot(sFilename_out)
                    logger.info(f'✓ Visualization saved to: {sFilename_out}')

                    # Verify file was created
                    if os.path.exists(sFilename_out):
                        file_size = os.path.getsize(sFilename_out)
                        logger.info(f'  File size: {file_size / 1024:.1f} KB')
                    else:
                        logger.warning(f'Screenshot command executed but file not found: {sFilename_out}')

                    plotter.close()
                    return True

                except Exception as e:
                    logger.error(f'Failed to save screenshot: {e}')
                    logger.error(f'Traceback: {traceback.format_exc()}')
                    plotter.close()
                    return False
            else:
                # Interactive display
                try:
                    logger.info('Opening interactive visualization window...')
                    plotter.show()
                    return True
                except Exception as e:
                    logger.error(f'Failed to display interactive visualization: {e}')
                    logger.error(f'Ensure display environment is available (X11, Wayland, etc.)')
                    logger.error(f'Traceback: {traceback.format_exc()}')
                    plotter.close()
                    return False

        except ImportError as e:
            logger.error(f'Missing required GeoVista dependencies: {e}')
            return False

        except Exception as e:
            logger.error(f'Unexpected error during target mesh visualization: {e}')
            logger.error(f'Error type: {type(e).__name__}')
            logger.error(f'Traceback: {traceback.format_exc()}')
            return False

    def _process_single_polygon(self, polygon, aFilename_source_raster, gdal_warp_options_base, feature_id):
        """
        Process a single polygon with GDAL Warp operation.

        Args:
            polygon: OGR Polygon geometry to clip raster
            aFilename_source_raster: List of source raster filenames
            gdal_warp_options_base: Base GDAL warp options dictionary
            feature_id: Feature identifier for logging

        Returns:
            tuple: (dataset, data_array, geotransform) or (None, None, None) on failure

        Raises:
            TimeoutError: If operation exceeds WARP_TIMEOUT_SECONDS
        """
        try:
            # Use WKT for faster performance
            pPolygonWKT = polygon.ExportToWkt()
            warp_options = gdal_warp_options_base.copy()
            warp_options['cutlineWKT'] = pPolygonWKT
            pWrapOption = gdal.WarpOptions(**warp_options)

            # Capture GDAL errors during Warp operation
            gdal.PushErrorHandler('CPLQuietErrorHandler')

            # Set a timeout for the operation (using alarm signal on Unix)
            def timeout_handler(signum, frame):
                raise TimeoutError(f"GDAL Warp operation timed out after {WARP_TIMEOUT_SECONDS} seconds for feature {feature_id}")

            if hasattr(signal, 'SIGALRM'):  # Unix systems only
                signal.signal(signal.SIGALRM, timeout_handler)
                signal.alarm(WARP_TIMEOUT_SECONDS)

            logger.debug(f"Starting GDAL Warp for single polygon feature {feature_id}")
            pDataset_clip_warped = gdal.Warp('', aFilename_source_raster, options=pWrapOption)

            if hasattr(signal, 'SIGALRM'):
                signal.alarm(0)  # Cancel the alarm

            gdal.PopErrorHandler()

            # Check if the operation succeeded
            if pDataset_clip_warped is None:
                gdal_error = gdal.GetLastErrorMsg()
                logger.error(f"GDAL Warp failed for single polygon feature {feature_id}: {gdal_error}")
                return None, None, None

            newGeoTransform = pDataset_clip_warped.GetGeoTransform()
            aData_clip = pDataset_clip_warped.ReadAsArray()

            # Check if data was successfully read
            if aData_clip is None:
                logger.error(f"Failed to read array data for single polygon feature {feature_id}")
                pDataset_clip_warped = None
                return None, None, None

            # Check for reasonable data dimensions
            if aData_clip.size == 0:
                logger.warning(f"Empty data array for single polygon feature {feature_id}")
                pDataset_clip_warped = None
                return None, None, None

            return pDataset_clip_warped, aData_clip, newGeoTransform

        except TimeoutError as e:
            logger.error(str(e))
            return None, None, None
        except Exception as e:
            logger.error(f"Unexpected exception during single polygon processing for feature {feature_id}: {str(e)}")
            return None, None, None

    def _process_multipolygon_idl(self, multipolygon, aFilename_source_raster, gdal_warp_options_base, feature_id, dMissing_value):
        """
        Process a multipolygon that crosses the International Date Line (IDL).

        Handles IDL-crossing features by processing each polygon part separately
        and merging the results into a single data array.

        Args:
            multipolygon: OGR MultiPolygon geometry that crosses IDL
            aFilename_source_raster: List of source raster filenames
            gdal_warp_options_base: Base GDAL warp options dictionary
            feature_id: Feature identifier for logging
            dMissing_value: Missing/NoData value for the raster

        Returns:
            tuple: (merged_data_array, merged_geotransform) or (None, None) on failure
            - merged_data_array: 1D array containing all valid data values
            - merged_geotransform: Geotransform from first polygon part
        """
        try:
            nGeometries = multipolygon.GetGeometryCount()
            logger.info(f"Processing {nGeometries} polygon parts for IDL-crossing feature {feature_id}")

            merged_datasets = []
            merged_data_arrays = []
            merged_transforms = []

            # Process each polygon part separately
            for iPart in range(nGeometries):
                polygon_part = multipolygon.GetGeometryRef(iPart)

                # Process this polygon part
                pDataset_part, aData_part, transform_part = self._process_single_polygon(
                    polygon_part, aFilename_source_raster, gdal_warp_options_base,
                    f"{feature_id}_part{iPart}")

                if pDataset_part is None:
                    logger.warning(f"Failed to process polygon part {iPart} of feature {feature_id}")
                    continue

                merged_datasets.append(pDataset_part)
                merged_data_arrays.append(aData_part)
                merged_transforms.append(transform_part)

            if not merged_data_arrays:
                logger.error(f"No polygon parts could be processed for IDL feature {feature_id}")
                return None, None

            # Merge the data arrays and transforms
            merged_data, merged_transform = self._merge_raster_parts(merged_data_arrays, merged_transforms, feature_id, dMissing_value)

            # Clean up datasets
            for dataset in merged_datasets:
                dataset = None

            return merged_data, merged_transform

        except Exception as e:
            logger.error(f"Error processing multipolygon IDL feature {feature_id}: {str(e)}")
            return None, None

    def _merge_raster_parts(self, data_arrays, transforms, feature_id, dMissing_value):
        """
        Merge multiple raster arrays from IDL-split polygons into a single data array.

        Extracts all valid (non-missing) data values from multiple raster arrays
        and concatenates them into a single 1D array for statistics calculation.

        Args:
            data_arrays: List of numpy arrays from different polygon parts
            transforms: List of geotransforms (one per data array)
            feature_id: Feature identifier for logging
            dMissing_value: Missing/NoData value to exclude from results

        Returns:
            tuple: (merged_1D_array, first_transform)
            - merged_1D_array: 1D array containing all valid data values
            - first_transform: Geotransform from first polygon part (dummy value)
        """
        try:
            if len(data_arrays) == 1:
                # Single array: extract valid data as 1D
                data_array = data_arrays[0]
                valid_mask = data_array != dMissing_value
                valid_data = data_array[valid_mask]
                return valid_data.flatten() if valid_data.size > 0 else np.array([dMissing_value]), transforms[0]

            # Multiple arrays: concatenate all valid data values
            all_valid_data = []

            for data_array in data_arrays:
                # Extract all valid (non-nodata) values
                valid_mask = data_array != dMissing_value
                valid_data = data_array[valid_mask]

                if valid_data.size > 0:
                    all_valid_data.append(valid_data.flatten())

            if not all_valid_data:
                logger.warning(f"No valid data found in any part for IDL feature {feature_id}")
                return np.array([dMissing_value]), transforms[0]

            # Concatenate all valid data into a single 1D array
            merged_data = np.concatenate(all_valid_data)

            logger.info(f"Successfully merged {len(data_arrays)} raster parts for feature {feature_id}: "
                       f"{merged_data.size} valid pixels")

            # Return 1D array of valid data only (no transform needed for statistics)
            return merged_data, transforms[0]

        except Exception as e:
            logger.error(f"Error merging raster parts for feature {feature_id}: {str(e)}")
            return None, None
