
# Extract module for uraster - contains remap workflow functions
import os
import logging
import traceback
import signal
import sys
import time
import threading
import gc
import numpy as np
from osgeo import gdal, ogr, osr
from multiprocessing import Pool, cpu_count
import multiprocessing
import tempfile
import shutil

from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
from pyearth.gis.geometry.international_date_line_utility import split_international_date_line_polygon_coordinates, check_cross_international_date_line_polygon
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename

# Try to import psutil for memory monitoring (optional)
try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False

from uraster.classes.sraster import sraster

# Set up logging
logger = logging.getLogger(__name__)
crs = "EPSG:4326"
pDriver_geojson = ogr.GetDriverByName('GeoJSON')
pDriver_shp = ogr.GetDriverByName('ESRI Shapefile')

# Constants for processing thresholds
IDL_LONGITUDE_THRESHOLD = 100  # Degrees - threshold for detecting IDL crossing
WARP_TIMEOUT_SECONDS = 30      # Seconds - timeout for GDAL Warp operations
PROGRESS_REPORT_INTERVAL = 5 # Report progress every N features
MAX_CONSECUTIVE_FAILURES = 10   # Maximum consecutive failures before stopping
HEARTBEAT_INTERVAL = 5          # Seconds between heartbeat logs during long operations

class TimeoutError(Exception):
    """Custom timeout exception for cross-platform compatibility"""
    pass

class CrossPlatformTimeout:
    """Cross-platform timeout mechanism using threading"""

    def __init__(self, timeout_seconds, operation_name="operation"):
        self.timeout_seconds = timeout_seconds
        self.operation_name = operation_name
        self.timer = None
        self.timed_out = False

    def __enter__(self):
        self.start_timer()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cancel_timer()
        if self.timed_out:
            raise TimeoutError(f"{self.operation_name} timed out after {self.timeout_seconds} seconds")

    def start_timer(self):
        """Start the timeout timer"""
        self.timer = threading.Timer(self.timeout_seconds, self._timeout_handler)
        self.timer.start()

    def cancel_timer(self):
        """Cancel the timeout timer"""
        if self.timer:
            self.timer.cancel()

    def _timeout_handler(self):
        """Handle timeout event"""
        self.timed_out = True
        logger.error(f"Timeout: {self.operation_name} exceeded {self.timeout_seconds} seconds")


def _warp_worker(src_files, warp_options_dict, tmp_filename):
    """
    Worker to run gdal.Warp in a child process and write output to tmp_filename.
    This must be a module-level function so multiprocessing can spawn it on Windows.
    Exits with code 0 on success, non-zero on failure.
    """
    try:
        # Import locally in child process to avoid pickling issues
        from osgeo import gdal
        import os
        opts = warp_options_dict.copy()
        # Ensure spatial reference entries are serializable (WKT strings)
        for key in ('dstSRS', 'srcSRS'):
            if key in opts and hasattr(opts[key], 'ExportToWkt'):
                try:
                    opts[key] = opts[key].ExportToWkt()
                except Exception:
                    opts[key] = str(opts[key])

        pWrapOption = gdal.WarpOptions(**opts)
        ds = gdal.Warp(tmp_filename, src_files, options=pWrapOption)
        if ds is None:
            try:
                if os.path.exists(tmp_filename):
                    os.remove(tmp_filename)
            except Exception:
                pass
            # non-zero exit shows failure
            sys.exit(2)
        # close dataset to flush to disk
        ds = None
        sys.exit(0)
    except Exception:
        try:
            if os.path.exists(tmp_filename):
                os.remove(tmp_filename)
        except Exception:
            pass
        sys.exit(3)


def _warp_with_timeout(src_files, warp_options_dict, timeout_seconds=WARP_TIMEOUT_SECONDS):
    """
    Run gdal.Warp in a separate process with a timeout. Returns (ndarray, geotransform) or (None, None).
    """
    # Ensure spatial references are serializable (WKT)
    opts = warp_options_dict.copy()
    for key in ('dstSRS', 'srcSRS'):
        if key in opts and hasattr(opts[key], 'ExportToWkt'):
            try:
                opts[key] = opts[key].ExportToWkt()
            except Exception:
                opts[key] = str(opts[key])

    tf = tempfile.NamedTemporaryFile(suffix='.tif', delete=False)
    tmp_filename = tf.name
    tf.close()

    proc = multiprocessing.Process(target=_warp_worker, args=(src_files, opts, tmp_filename))
    proc.start()
    proc.join(timeout_seconds)

    if proc.is_alive():
        try:
            proc.terminate()
            proc.join(5)
        except Exception:
            pass
        try:
            os.remove(tmp_filename)
        except Exception:
            pass
        logger.error(f"GDAL Warp timed out after {timeout_seconds}s and was terminated")
        return None, None

    if proc.exitcode != 0:
        try:
            if os.path.exists(tmp_filename):
                os.remove(tmp_filename)
        except Exception:
            pass
        logger.error(f"GDAL Warp process failed with exit code {proc.exitcode}")
        return None, None

    # Read the result produced by the worker
    try:
        ds = gdal.Open(tmp_filename, gdal.GA_ReadOnly)
        if ds is None:
            logger.error("Child gdal.Warp created no readable output")
            try:
                os.remove(tmp_filename)
            except Exception:
                pass
            return None, None
        newGeoTransform = ds.GetGeoTransform()
        aData_clip = ds.ReadAsArray()
        ds = None
        try:
            os.remove(tmp_filename)
        except Exception:
            pass
        if aData_clip is None:
            logger.error("Failed to read warped data from temp file")
            return None, None
        return aData_clip, newGeoTransform
    except Exception as e:
        logger.error(f"Error reading warp output: {e}")
        try:
            os.remove(tmp_filename)
        except Exception:
            pass
        return None, None


def _process_single_polygon(polygon, aFilename_source_raster, gdal_warp_options_base, feature_id, iFlag_verbose=False):
    """
    Process a single polygon with GDAL Warp operation.

    Args:
        polygon: OGR Polygon geometry to clip raster
        aFilename_source_raster: List of source raster filenames
        gdal_warp_options_base: Base GDAL warp options dictionary
        feature_id: Feature identifier for logging
        iFlag_verbose (bool, optional): If True, print detailed progress messages.
            If False, only print error messages. Default is False.

    Returns:
        tuple: (dataset, data_array, geotransform) or (None, None, None) on failure

    Raises:
        TimeoutError: If operation exceeds WARP_TIMEOUT_SECONDS
    """
    srs_wgs84 = None
    pDataset_clip = None
    pLayer_clip = None
    pFeature_clip = None
    pPolygonWKT_file = None
    try:
        # Use WKT for faster performance
        srs_wgs84 = osr.SpatialReference()
        srs_wgs84.ImportFromEPSG(4326)
        pPolygonWKT = polygon.ExportToWkt()
        #make a copy of the warp options to modify
        gdal_warp_options = gdal_warp_options_base.copy()
        #save WKT to a temporary file using the memory driver

        pPolygonWKT_file = '/vsimem/polygon_wkt_' + str(feature_id) + '.shp'
        pDataset_clip = pDriver_shp.CreateDataSource(pPolygonWKT_file)
        pLayer_clip = pDataset_clip.CreateLayer('polygon_layer', geom_type=ogr.wkbPolygon, srs=srs_wgs84)
        pFeature_clip = ogr.Feature(pLayer_clip.GetLayerDefn())
        pFeature_clip.SetGeometry(polygon)
        pLayer_clip.CreateFeature(pFeature_clip)
        pDataset_clip.FlushCache()

        gdal_warp_options['cutlineDSName'] = pPolygonWKT_file

        # Run GDAL Warp in a separate process with enforced timeout to avoid blocking the main process
        warp_start_time = time.time()
        aData_clip, newGeoTransform = _warp_with_timeout(aFilename_source_raster, gdal_warp_options, timeout_seconds=WARP_TIMEOUT_SECONDS)
        warp_duration = time.time() - warp_start_time

        if aData_clip is None:
            logger.error(f"GDAL Warp failed or timed out for single polygon feature {feature_id}")
            return None, None, None

        # Check if data was successfully read
        if aData_clip is None:
            logger.error(f"Failed to read array data for single polygon feature {feature_id}")
            return None, None, None

        # Check for reasonable data dimensions
        if aData_clip.size == 0:
            logger.warning(f"Empty data array for single polygon feature {feature_id}")
            return None, None, None

        # Make a copy of the data to allow immediate dataset cleanup
        aData_clip_copy = aData_clip.copy()

        # No dataset object in parent process; data array is returned from worker

        return aData_clip_copy, newGeoTransform  # Return None for dataset since it's cleaned up

    except TimeoutError as e:
        logger.error(str(e))
        return None, None
    except Exception as e:
        logger.error(f"Unexpected exception during single polygon processing for feature {feature_id}: {str(e)}")
        return None, None
    finally:
        # Clean up all spatial reference and OGR objects to prevent memory leaks
        if pPolygonWKT_file is not None:
            try:
                gdal.Unlink(pPolygonWKT_file)
            except Exception:
                pass
        if srs_wgs84 is not None:
            srs_wgs84 = None
        if pFeature_clip is not None:
            pFeature_clip = None
        if pLayer_clip is not None:
            pLayer_clip = None
        if pDataset_clip is not None:
            pDataset_clip = None


def _process_multipolygon_idl(multipolygon, aFilename_source_raster, gdal_warp_options_base, feature_id, dMissing_value, iFlag_verbose=False):
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
        iFlag_verbose (bool, optional): If True, print detailed progress messages.
            If False, only print error messages. Default is False.

    Returns:
        tuple: (merged_data_array, merged_geotransform) or (None, None) on failure
        - merged_data_array: 1D array containing all valid data values
        - merged_geotransform: Geotransform from first polygon part
    """
    try:
        nGeometries = multipolygon.GetGeometryCount()
        if iFlag_verbose:
            logger.info(f"Processing {nGeometries} polygon parts for IDL-crossing feature {feature_id}")


        merged_data_arrays = []
        merged_transforms = []

        # Process each polygon part separately
        for iPart in range(nGeometries):
            polygon_part = multipolygon.GetGeometryRef(iPart)

            # Process this polygon part
            aData_part, transform_part = _process_single_polygon(
                polygon_part, aFilename_source_raster, gdal_warp_options_base,
                f"{feature_id}_part{iPart}", iFlag_verbose)

            if aData_part is None:
                logger.warning(f"Failed to process polygon part {iPart} of feature {feature_id}")
                continue

            merged_data_arrays.append(aData_part)
            merged_transforms.append(transform_part)

        if not merged_data_arrays:
            logger.error(f"No polygon parts could be processed for IDL feature {feature_id}")
            return None, None

        # Merge the data arrays and transforms
        merged_data, merged_transform = _merge_raster_parts(merged_data_arrays, merged_transforms, feature_id, dMissing_value, iFlag_verbose)

        return merged_data, merged_transform

    except Exception as e:
        logger.error(f"Error processing multipolygon IDL feature {feature_id}: {str(e)}")
        return None, None


def _merge_raster_parts(data_arrays, transforms, feature_id, dMissing_value, iFlag_verbose=False):
    """
    Merge multiple raster arrays from IDL-split polygons into a single data array.

    Extracts all valid (non-missing) data values from multiple raster arrays
    and concatenates them into a single 1D array for statistics calculation.

    Args:
        data_arrays: List of numpy arrays from different polygon parts
        transforms: List of geotransforms (one per data array)
        feature_id: Feature identifier for logging
        dMissing_value: Missing/NoData value to exclude from results
        iFlag_verbose (bool, optional): If True, print detailed progress messages.
            If False, only print error messages. Default is False.

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

        if iFlag_verbose:
            logger.info(f"Successfully merged {len(data_arrays)} raster parts for feature {feature_id}: "
                       f"{merged_data.size} valid pixels")

        # Return 1D array of valid data only (no transform needed for statistics)
        return merged_data, transforms[0]

    except Exception as e:
        logger.error(f"Error merging raster parts for feature {feature_id}: {str(e)}")
        return None, None


def _process_feature_with_enhanced_handling(self, pPolygon, aFilename_source_raster, gdal_warp_options_base, feature_id, dMissing_value, iFlag_verbose=False):
    """
    Process a feature with enhanced error handling, timeout, and memory management.

    Args:
        self: uraster instance for accessing methods
        pPolygon: OGR Polygon or MultiPolygon geometry
        aFilename_source_raster: List of source raster filenames
        gdal_warp_options_base: Base GDAL warp options dictionary
        feature_id: Feature identifier for logging
        dMissing_value: Missing/NoData value for the raster
        iFlag_verbose (bool, optional): If True, print detailed progress messages.

    Returns:
        tuple: (dataset, data_array, geotransform) or (None, None, None) on failure
    """
    try:
        # Check geometry type and process accordingly
        sGeometry_type = pPolygon.GetGeometryName()

        if sGeometry_type == "POLYGON":
            # Use cross-platform timeout for single polygon processing
            with CrossPlatformTimeout(WARP_TIMEOUT_SECONDS, f"GDAL Warp for feature {feature_id}"):
                return _process_single_polygon(
                    pPolygon, aFilename_source_raster, gdal_warp_options_base, feature_id, iFlag_verbose)

        elif sGeometry_type == "MULTIPOLYGON":
            # Process multipolygon with timeout
            with CrossPlatformTimeout(WARP_TIMEOUT_SECONDS * 2, f"GDAL Warp for multipolygon feature {feature_id}"):
                merged_data, merged_transform = _process_multipolygon_idl(
                    pPolygon, aFilename_source_raster, gdal_warp_options_base, feature_id, dMissing_value, iFlag_verbose)

                if merged_data is None:
                    return None, None, None

                return True, merged_data, merged_transform  # Return flag to indicate multipolygon processing
        else:
            logger.error(f"Unsupported geometry type for feature {feature_id}: {sGeometry_type}")
            return None, None, None

    except TimeoutError as e:
        logger.error(f"Timeout processing feature {feature_id}: {e}")
        return None, None, None
    except Exception as e:
        logger.error(f"Error processing feature {feature_id}: {e}")
        return None, None, None

def get_polygon_list(sFilename_source_mesh, iFlag_verbose=False):
    aPolygon = []

    # Get the spatial reference of the mesh vector file
    pDataset_mesh = ogr.Open(sFilename_source_mesh, 0 ) #0 means read-only. 1 means writeable.
    if pDataset_mesh is None:
        logger.error(f"Failed to open mesh dataset: {sFilename_source_mesh}")
        return
    pLayer_mesh = pDataset_mesh.GetLayer(0)
    if pLayer_mesh is None:
        logger.error("Failed to get layer from mesh dataset")
        return
    nFeature = pLayer_mesh.GetFeatureCount()
    if iFlag_verbose:
        logger.info(f"Found {nFeature} features in mesh dataset")
    pSpatialRef_target = pLayer_mesh.GetSpatialRef()

    pSpatialRef_target_wkt = pSpatialRef_target.ExportToWkt() if pSpatialRef_target else None
    pLayer_mesh.ResetReading()
    pFeature_mesh = pLayer_mesh.GetNextFeature()
    while pFeature_mesh is not None:
        # Check if feature crosses the international date line
        pPolygon = pFeature_mesh.GetGeometryRef()
        sGeometry_type = pPolygon.GetGeometryName()
        # Read cellid from current feature
        current_cellid = pFeature_mesh.GetField('cellid')
        if sGeometry_type == "POLYGON":
            aCoord = get_geometry_coordinates(pPolygon)
            #first check whether a geometry crosses the IDL
            if check_cross_international_date_line_polygon(aCoord):  # Use constant
                # Create multipolygon to handle IDL crossing
                if iFlag_verbose:
                    logger.info(f'Feature {i} crosses the international date line, splitting into multiple parts.')
                pMultipolygon = ogr.Geometry(ogr.wkbMultiPolygon)
                aCoord_gcs_split = split_international_date_line_polygon_coordinates(aCoord)

                for aCoord_gcs in aCoord_gcs_split:
                    #create a polygon (not just ring) and add it to the multipolygon
                    ring = ogr.Geometry(ogr.wkbLinearRing)
                    for iCoord in range(aCoord_gcs.shape[0]):
                        ring.AddPoint(aCoord_gcs[iCoord, 0], aCoord_gcs[iCoord, 1])

                    ring.CloseRings()
                    # Create polygon from ring
                    polygon_part = ogr.Geometry(ogr.wkbPolygon)
                    polygon_part.AddGeometry(ring)
                    # check validity
                    if not polygon_part.IsValid():
                        polygon_part.FlattenTo2D()
                        print(polygon_part.ExportToWkt())
                        logger.warning(f'Polygon part of feature {i} is invalid after splitting at IDL.')

                    pMultipolygon.AddGeometry(polygon_part)  # Add polygon, not ring

                #create a feature from the multipolygon and add it to the list
                wkt = pMultipolygon.ExportToWkt()
                aPolygon.append((current_cellid, wkt))
                pFeature_mesh = pLayer_mesh.GetNextFeature()
            else:
                wkt = pPolygon.ExportToWkt()
                aPolygon.append((current_cellid, wkt))
                pFeature_mesh = pLayer_mesh.GetNextFeature()
        else:
            #a feature may still have multipolygon
            wkt = pPolygon.ExportToWkt()
            aPolygon.append((current_cellid, wkt))
            pFeature_mesh = pLayer_mesh.GetNextFeature()

        i += 1

        # Progress reporting during feature pre-processing
        if i % 1000 == 0:
            logger.info(f"Pre-processed {i} features...")

    if iFlag_verbose:
        logger.info(f"run_remap: Pre-processing completed. Found {len(aPolygon)} features to process")

    return aPolygon

def run_remap(sFilename_target_mesh,
                sFilename_source_mesh,
                aFilename_source_raster,
                iFlag_remap_method_in = 1,
              iFlag_stat_in = 1,
              iFlag_save_clipped_raster_in=0,
              sFolder_raster_out_in = None,
              sFormat_in='GTiff',
              iFlag_verbose=False):
    """
    Perform zonal statistics by clipping raster data to mesh polygons.

    Main processing method that extracts raster values for each mesh cell polygon
    and computes statistics (mean, min, max, std, sum, count).

    Args:

        sFilename_vector_out (str): Output vector file path with computed statistics
        sFilename_source_mesh_in (str, optional): Input mesh polygon file.
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
        iFlag_verbose (bool, optional): If True, print detailed progress messages.
            If False, only print error messages. Default is False.

    Returns:
        None

    Note:
        - Handles IDL-crossing polygons automatically
        - Generates failure report for problematic features
        - Supports multiple input rasters (uses first in list)
    """
    if iFlag_remap_method_in not in [1, 2, 3]:
        logger.error("Invalid remap method specified. Must be 1 (nearest), 2 (bilinear), or 3 (weighted average).")
        return
    iFlag_remap_method = iFlag_remap_method_in


    if iFlag_verbose:
        logger.info("run_remap: Starting input file validation...")
    #check input files
    for idx, sFilename_raster in enumerate(aFilename_source_raster):
        if iFlag_verbose:
            logger.info(f"Checking raster file {idx+1}/{len(aFilename_source_raster)}: {os.path.basename(sFilename_raster)}")
        if os.path.exists(sFilename_raster):
            pass
        else:
            logger.error('The raster file does not exist!')
            return

    if iFlag_verbose:
        logger.info(f"Checking source mesh file: {os.path.basename(sFilename_source_mesh)}")
    if os.path.exists(sFilename_source_mesh):
        pass
    else:
        logger.error('The vector mesh file does not exist!')
        return
    if iFlag_verbose:
        logger.info("Input file validation completed successfully")

    # Determine output vector format from filename extension
    pDriver_vector = get_vector_driver_from_filename(sFilename_target_mesh)
    #check the input raster data format and decide gdal driver
    if sFormat_in is not None:
        sDriverName = sFormat_in
    else:
        sDriverName = 'GTiff'

    if os.path.exists(sFilename_target_mesh):
        #remove the file using the vector driver
        pDriver_vector.DeleteDataSource(sFilename_target_mesh)

    pDriver_raster = gdal.GetDriverByName(sDriverName)

    #get the raster file extension
    sFilename_raster = aFilename_source_raster[0]  #just use the first raster to get the extension
    sExtension = os.path.splitext(sFilename_raster)[1]
    sName = os.path.basename(sFilename_raster)
    sRasterName_no_extension = os.path.splitext(sName)[0]

    if iFlag_verbose:
        logger.info("run_remap: Reading raster metadata and determining processing bounds...")
    #use the sraster class the check the raster
    dX_left = -180.0
    dX_right = 180.0
    dY_top = 90.0
    dY_bot = -90.0
    dPixelWidth = None
    pPixelHeight = None
    #narrow the range to speed up the processing
    #also get the highest resolution
    for idx, sFilename_raster in enumerate(aFilename_source_raster):
        if iFlag_verbose:
            logger.info(f"Reading metadata for raster {idx+1}/{len(aFilename_source_raster)}: {os.path.basename(sFilename_raster)}")
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
    sRemap_method_auto, iRemap_method_auto = _determine_optimal_resampling(dPixelWidth, abs(pPixelHeight), iFlag_verbose)

    # Use automatically determined method if it's more conservative than user setting
    # Priority: weighted averaging > nearest neighbor
    if iRemap_method_auto == 3 and iFlag_remap_method != 3:
        logger.warning(f"Overriding user remap method ({iFlag_remap_method}) with automatic selection (3 - weighted average)")
        logger.warning("This is necessary due to mesh/raster resolution compatibility")
        sRemap_method = sRemap_method_auto
    else:
        # Use user's preferred method
        if iFlag_remap_method == 1:
            sRemap_method = 'near'
        elif iFlag_remap_method == 2:
            sRemap_method = 'near'
        elif iFlag_remap_method == 3:
            sRemap_method = 'average'
        if iFlag_verbose:
            logger.info(f"Using user-specified remap method: {sRemap_method}")


    if iFlag_verbose:
        logger.info("run_remap: Opening mesh dataset and analyzing features...")


    aPolygon = get_polygon_list(sFilename_source_mesh, iFlag_verbose)

    options = ['COMPRESS=DEFLATE', 'PREDICTOR=2']
    sDriverName = 'MEM'

    #create a polygon feature to save the output
    pDataset_out = pDriver_vector.CreateDataSource(sFilename_target_mesh)
    pLayer_out = pDataset_out.CreateLayer('uraster', pSpatialRef_target, ogr.wkbPolygon)
    pLayer_defn_out = pLayer_out.GetLayerDefn()
    pFeature_out = ogr.Feature(pLayer_defn_out)

    #add id, area and mean, min, max, std of the raster
    pLayer_out.CreateField(ogr.FieldDefn('cellid', ogr.OFTInteger))
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

    i = 1
    consecutive_failures = 0

    if iFlag_verbose:
        logger.info("run_remap: Pre-fetching features and analyzing geometries...")



    start_time = time.time()

    # Add crash tracking variables
    failed_features = []
    successful_features = 0
    logger.info("run_remap: Starting main feature processing loop...")


    #reset i
    i = 1

    # Process features with enhanced error handling and memory management
    for idx, pFeature_mesh in enumerate(aFeatures):
        try:
            # Enhanced progress reporting for debugging (reduce frequency for large datasets)
            if i == 1 or i % PROGRESS_REPORT_INTERVAL == 0:
                elapsed = time.time() - start_time
                rate = i / elapsed if elapsed > 0 else 0
                if iFlag_verbose:
                    logger.info(f"Starting feature {i}/{len(aFeatures)} (Rate: {rate:.2f} features/sec)")

            logger.debug(f"Processing feature {i}: Getting geometry and envelope...")
            sClip = f"{i:08d}"  # f-string is faster than format
            pPolygon = pFeature_mesh.GetGeometryRef()
            if pPolygon is None:
                logger.error(f"Feature {i} has no geometry!")
                consecutive_failures += 1
                i += 1
                continue

            # Fast envelope check first
            logger.debug(f"Feature {i}: Getting envelope...")
            minX, maxX, minY, maxY = pPolygon.GetEnvelope()
            logger.debug(f"Feature {i}: Envelope = ({minX:.3f}, {maxX:.3f}, {minY:.3f}, {maxY:.3f})")
            if (minX > dX_right or maxX < dX_left or
                minY > dY_top or maxY < dY_bot or
                not pPolygon or pPolygon.IsEmpty() or not pPolygon.IsValid()):
                i += 1
                consecutive_failures += 1
                print(f"Skipping feature {i} due to envelope check.")
                print(pPolygon.ExportToWkt())
                continue

            # Process valid polygons with circuit breaker protection
            feature_start_time = time.time()

            try:
                # Process feature directly without circuit breaker
                result = _process_feature_with_enhanced_handling(
                    pPolygon, aFilename_source_raster, gdal_warp_options_base,
                    i, dMissing_value, iFlag_verbose)
                aData_clip, newGeoTransform = result

                # Check if processing failed - aData_clip is None indicates failure
                if aData_clip is None:
                    error_msg = f"Failed to process feature {i}"
                    logger.error(error_msg)
                    failed_features.append({'feature_id': i, 'error': error_msg, 'envelope': (minX, maxX, minY, maxY)})
                    consecutive_failures += 1
                    i += 1
                    continue

                # Reset consecutive failures on success
                consecutive_failures = 0
                successful_features += 1

            except Exception as processing_error:
                error_msg = f"Processing error for feature {i}: {str(processing_error)}"
                logger.error(error_msg)
                failed_features.append({'feature_id': i, 'error': error_msg, 'envelope': (minX, maxX, minY, maxY)})
                consecutive_failures += 1
                i += 1
                continue

            # Calculate area early - this is independent of raster processing
            aCoords_gcs = get_geometry_coordinates(pPolygon)
            sGeometry_type = pPolygon.GetGeometryName()
            if sGeometry_type == 'POLYGON':
                dArea = calculate_polygon_area(aCoords_gcs[:,0], aCoords_gcs[:,1])
            else:
                dArea = 0.0
                for iPart in range(pPolygon.GetGeometryCount()):
                    pPolygon_part = pPolygon.GetGeometryRef(iPart)
                    aCoords_part = get_geometry_coordinates(pPolygon_part)
                    dArea += calculate_polygon_area(aCoords_part[:,0], aCoords_part[:,1])

            logger.debug(f"Successfully processed feature {i} in {time.time() - feature_start_time:.2f}s")

        except Exception as e:
            error_msg = f"Unexpected exception during feature {i} processing: {str(e)}"
            logger.error(error_msg)
            logger.error(f"Exception type: {type(e).__name__}")
            logger.error(f"Traceback: {traceback.format_exc()}")
            failed_features.append({'feature_id': i, 'error': error_msg, 'envelope': (minX, maxX, minY, maxY)})
            consecutive_failures += 1
            i += 1
            continue

        # Check if too many consecutive failures occurred
        if consecutive_failures >= MAX_CONSECUTIVE_FAILURES:
            logger.error(f"Too many consecutive failures ({consecutive_failures}). Stopping processing.")
            break

        # Optimize data type conversion and missing value handling
        if aData_clip is not None:
            # Use numpy's faster operations
            aData_clip = aData_clip.astype(np.int32, copy=False)  # Avoid unnecessary copy
            np.place(aData_clip, aData_clip == dMissing_value, -9999)  # Faster than boolean indexing

            if iFlag_save_clipped_raster_in == 1 :
                if sGeometry_type == ogr.wkbPolygon : #we cannot save the multipolygon directly
                    sFilename_raster_out = os.path.join(sFolder_raster_out_in, f'{sRasterName_no_extension}_clip_{sClip}{sExtension}')
                    # Delete existing file if it exists (GDAL Create() doesn't overwrite)
                    if os.path.exists(sFilename_raster_out):
                        try:
                            pDriver_raster.Delete(sFilename_raster_out)
                        except:
                            os.remove(sFilename_raster_out)  # Fallback if GDAL delete fails
                    iNewWidth = aData_clip.shape[1]
                    iNewHeigh = aData_clip.shape[0]
                    pDataset_clip = pDriver_raster.Create(sFilename_raster_out, iNewWidth, iNewHeigh, 1, eType , options= options)
                    pDataset_clip.SetGeoTransform( newGeoTransform )
                    pDataset_clip.SetProjection( pSpatialRef_target_wkt)
                    #set the no data value
                    pBand = pDataset_clip.GetRasterBand(1)
                    pBand.SetNoDataValue(-9999)
                    pBand.WriteArray(aData_clip)
                    pBand.FlushCache()  # Corrected method name to FlushCache()
                    pDataset_clip.FlushCache()
                    pDataset_clip = None

            #create a polygon feature to save the output
            pFeature_out.SetGeometry(pPolygon.Clone())
            # Use the pre-read cellid from the parallel array
            actual_cellid = aCellIDs_features[idx]
            pFeature_out.SetField('cellid', actual_cellid)
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

            # Explicit cleanup - no dataset cleanup needed since _process_single_polygon
            # handles immediate cleanup for memory management
            pass

            # Clean up data array after use
            if 'aData_clip' in locals():
                del aData_clip

            # Progress reporting and memory management
            if i % PROGRESS_REPORT_INTERVAL == 0:
                elapsed = time.time() - start_time
                rate = i / elapsed
                eta = (len(aFeatures) - i) / rate if rate > 0 else 0
                #if iFlag_verbose:
                print(f"Processed {i}/{len(aFeatures)} features ({rate:.2f} features/sec, ETA: {eta:.0f}s)")

                # Force garbage collection every PROGRESS_REPORT_INTERVAL features
                gc.collect()

                # Optional: Log memory usage if psutil is available
                if PSUTIL_AVAILABLE and iFlag_verbose:
                    process = psutil.Process()
                    memory_mb = process.memory_info().rss / 1024 / 1024
                    logger.info(f"Memory usage: {memory_mb:.1f} MB")

        i += 1


    pDataset_out.FlushCache()  # Flush once after all features are added
    pDataset_out = None        # Close the dataset
    pDataset_mesh = None

    # Clean up spatial reference objects to prevent memory leaks
    pSpatialRef_target = None

    # Report processing summary
    total_time = time.time() - start_time
    if iFlag_verbose:
        logger.info(f"Processing completed in {total_time:.2f} seconds")
        logger.info(f"Successfully processed: {successful_features} features")
        logger.info(f"Failed features: {len(failed_features)}")

    if failed_features:
        if iFlag_verbose:
            logger.warning("Failed features summary:")
            for failed in failed_features[:10]:  # Show first 10 failures
                logger.warning(f"  Feature {failed['feature_id']}: {failed['error']}")
            if len(failed_features) > 10:
                logger.warning(f"  ... and {len(failed_features) - 10} more failures")

        # Save failure report to file
        # Generate failure report filename by replacing extension with '_failures.log'
        base_name = os.path.splitext(sFilename_target_mesh_out)[0]
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
            if iFlag_verbose:
                logger.info(f"Failure report saved to: {failure_report_file}")
        except Exception as e:
            logger.error(f"Could not save failure report: {e}")

    return

def _determine_optimal_resampling(self, dPixelWidth, dPixelHeight, iFlag_verbose=False):
    """
    Determine optimal resampling method based on mesh and raster resolution comparison.
    Compares the characteristic mesh cell size with raster resolution to decide
    whether to use nearest neighbor (when raster is much finer) or weighted
    averaging (when mesh and raster resolutions are comparable).
    Args:
        dPixelWidth (float): Raster pixel width in degrees
        dPixelHeight (float): Raster pixel height in degrees (absolute value)
        iFlag_verbose (bool, optional): If True, print detailed progress messages.
            If False, only print error messages. Default is False.
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
    if iFlag_verbose:
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
        if iFlag_verbose:
            logger.warning(f"Mesh resolution is close to raster resolution (ratio: {dResolution_ratio:.2f})")
            logger.warning(f"Switching to WEIGHTED AVERAGING (average) for accuracy")
            logger.warning("Consider using higher resolution raster data for better results")
    else:
        # Raster is much finer than mesh - nearest neighbor is appropriate
        recommended_method = 'near'
        recommended_code = 1
        if iFlag_verbose:
            logger.info(f"Raster is significantly finer than mesh (ratio: {dResolution_ratio:.2f})")
            logger.info(f"Using NEAREST NEIGHBOR resampling (sufficient for this resolution ratio)")
    if iFlag_verbose:
        logger.info("="*60)
    return recommended_method, recommended_code