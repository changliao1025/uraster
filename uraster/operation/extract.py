
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
from concurrent.futures import ProcessPoolExecutor, as_completed
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

def _determine_optimal_resampling(dArea_mean, dPixelWidth, dPixelHeight, iFlag_verbose=False, dResolution_ratio_threshold=3.0):
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

    # Estimate characteristic mesh cell dimension (approximate square root of mean area)
    dMesh_characteristic_size = np.sqrt(dArea_mean)
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
        logger.info(f"  Threshold for weighted averaging: {dResolution_ratio_threshold:.2f}")
    # Decision logic
    if dResolution_ratio < dResolution_ratio_threshold:
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


def _process_single_polygon(feature_idx, wkt, aFilename_source_raster, gdal_warp_options_base, cellid, dMissing_value, iFlag_verbose=False):
    """
    Process a single polygon with GDAL Warp operation.

    Args:
        wkt: Well-Known Text representation of the polygon geometry to clip raster
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
    pPolygonWKT_file = None
    srs_wgs84 = None
    pFeature_clip = None
    pLayer_clip = None
    pDataset_clip = None
    try:
        # Use WKT for faster performance
        srs_wgs84 = osr.SpatialReference()
        srs_wgs84.ImportFromEPSG(4326)
        polygon = ogr.CreateGeometryFromWkt(wkt)
        #make a copy of the warp options to modify
        gdal_warp_options = gdal_warp_options_base.copy()
        #save WKT to a temporary file using the memory driver
        pPolygonWKT_file = '/vsimem/polygon_wkt_' + str(cellid) + '.shp'
        pDataset_clip = pDriver_shp.CreateDataSource(pPolygonWKT_file)
        pLayer_clip = pDataset_clip.CreateLayer('polygon', geom_type=ogr.wkbPolygon, srs=srs_wgs84)
        pFeature_clip = ogr.Feature(pLayer_clip.GetLayerDefn())
        pFeature_clip.SetGeometry(polygon)
        pLayer_clip.CreateFeature(pFeature_clip)
        pDataset_clip.FlushCache()

        gdal_warp_options['cutlineDSName'] = pPolygonWKT_file

        # Run GDAL Warp in a separate process with enforced timeout to avoid blocking the main process
        warp_start_time = time.time()
        pWrapOption = gdal.WarpOptions(**gdal_warp_options)
        pDataset_warp = gdal.Warp("", aFilename_source_raster, options=pWrapOption)
        aData_clip = pDataset_warp.ReadAsArray()
        newGeoTransform = pDataset_warp.GetGeoTransform()
        warp_duration = time.time() - warp_start_time
        if iFlag_verbose:
            logger.info(f"GDAL Warp completed in {warp_duration:.2f} seconds")

        if aData_clip is None:
            error_msg = f"GDAL Warp failed or timed out for single polygon feature {cellid}"
            logger.error(error_msg)
            return feature_idx, cellid, False, error_msg

        # Check for reasonable data dimensions
        if aData_clip.size == 0:
            error_msg = f"Empty data array for single polygon feature {cellid}"
            logger.warning(error_msg)
            return feature_idx, cellid, False, error_msg

        # Make a copy of the data to allow immediate dataset cleanup
        aData_clip_copy = aData_clip.copy()

        # Calculate statistics
        # Filter out missing/nodata values
        valid_data = aData_clip_copy[aData_clip_copy != dMissing_value]
        if len(valid_data) == 0:
            # No valid pixels for this feature: not treated as an error.
            # Return NaN stats and count 0 so caller can handle missing-data features.
            stats = {
                'mean': float(np.nan),
                'min': float(np.nan),
                'max': float(np.nan),
                'std': float(np.nan),
                'count': 0
            }
        else:
            # Compute statistics on valid data
            stats = {
                'mean': float(np.mean(valid_data)),
                'min': float(np.min(valid_data)),
                'max': float(np.max(valid_data)),
                'std': float(np.std(valid_data)),
                'count': int(len(valid_data))
            }
    finally:
        srs_wgs84 = None

    return feature_idx, cellid, True, stats

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

def _process_multipolygon_idl(feature_idx, wkt, aFilename_source_raster, gdal_warp_options_base, feature_id, dMissing_value, iFlag_verbose=False):
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
        multipolygon = ogr.CreateGeometryFromWkt(wkt)
        nGeometries = multipolygon.GetGeometryCount()
        if iFlag_verbose:
            logger.info(f"Processing {nGeometries} polygon parts for IDL-crossing feature {feature_id}")


        merged_data_arrays = []
        merged_transforms = []

        # Process each polygon part separately
        for iPart in range(nGeometries):
            polygon_part = multipolygon.GetGeometryRef(iPart)

            # Process this polygon part - convert geometry to WKT
            polygon_part_wkt = polygon_part.ExportToWkt()
            part_result = _process_single_polygon(
                feature_idx, polygon_part_wkt, aFilename_source_raster, gdal_warp_options_base,
                f"{feature_id}_part{iPart}", dMissing_value, iFlag_verbose)

            # _process_single_polygon returns (feature_idx, cellid, success, stats_or_error)
            if len(part_result) != 4 or not part_result[2]:
                logger.warning(f"Failed to process polygon part {iPart} of feature {feature_id}")
                continue

            # Extract the stats from the result and convert to data array for merging
            part_stats = part_result[3]
            if isinstance(part_stats, dict) and 'count' in part_stats:
                # Create a dummy array with the mean value repeated 'count' times
                # This preserves the statistical properties for merging
                part_data = np.full(part_stats['count'], part_stats['mean'])
                merged_data_arrays.append(part_data)
                merged_transforms.append(None)  # Transform not needed for statistics
            else:
                logger.warning(f"Invalid stats returned for polygon part {iPart} of feature {feature_id}")
                continue

        if not merged_data_arrays:
            logger.error(f"No polygon parts could be processed for IDL feature {feature_id}")
            return None, None

        # Merge the data arrays and transforms
        merged_data, merged_transform = _merge_raster_parts(merged_data_arrays, merged_transforms, feature_id, dMissing_value, iFlag_verbose)

        return merged_data, merged_transform

    except Exception as e:
        logger.error(f"Error processing multipolygon IDL feature {feature_id}: {str(e)}")
        return None, None

def _process_task(args):
    """
    Module-level worker for multiprocessing. Accepts a tuple:
    (feature_idx, cellid, wkt, aFilename_source_raster, gdal_warp_options_serial, dMissing_value, iFlag_verbose)
    Returns a tuple (feature_idx, cellid, success, stats_dict_or_error)
    """
    feature_idx, cellid, wkt, aFilename_source_raster, gdal_warp_options_base, dMissing_value, iFlag_verbose = args
    try:
        # Check geometry type and process accordingly
        pPolygon = ogr.CreateGeometryFromWkt(wkt)
        sGeometry_type = pPolygon.GetGeometryName()

        if sGeometry_type == "POLYGON":
            return _process_single_polygon(feature_idx,
                    wkt, aFilename_source_raster, gdal_warp_options_base, cellid, dMissing_value, iFlag_verbose)

        elif sGeometry_type == "MULTIPOLYGON":
            merged_data, merged_transform = _process_multipolygon_idl(feature_idx,
                    wkt, aFilename_source_raster, gdal_warp_options_base, cellid, dMissing_value, iFlag_verbose)

            if merged_data is None:
                return feature_idx, cellid, False, "Failed to process multipolygon"

            # Calculate statistics for multipolygon data
            valid_data = merged_data[merged_data != dMissing_value]
            if len(valid_data) == 0:
                return feature_idx, cellid, False, "No valid data found in multipolygon"

            stats = {
                'mean': float(np.mean(valid_data)),
                'min': float(np.min(valid_data)),
                'max': float(np.max(valid_data)),
                'std': float(np.std(valid_data)),
                'count': len(valid_data)
            }
            return feature_idx, cellid, True, stats
        else:
            logger.error(f"Unsupported geometry type for feature {cellid}: {sGeometry_type}")
            return feature_idx, cellid, False, f"Unsupported geometry type: {sGeometry_type}"

    except TimeoutError as e:
        logger.error(f"Timeout processing feature {cellid}: {e}")
        return feature_idx, cellid, False, f"Timeout: {str(e)}"
    except Exception as e:
        logger.error(f"Error processing feature {cellid}: {e}")
        return feature_idx, cellid, False, f"Error: {str(e)}"

def get_polygon_list(sFilename_source_mesh, iFlag_verbose=False):
    if iFlag_verbose:
        logger.info("run_remap: Pre-fetching features and analyzing geometries...")
    aPolygon = []
    aArea = []

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
    pSpatialRef_source = pLayer_mesh.GetSpatialRef()

    sProjection_source_wkt = pSpatialRef_source.ExportToWkt() if pSpatialRef_source else None
    pLayer_mesh.ResetReading()
    pFeature_mesh = pLayer_mesh.GetNextFeature()
    #reset i
    i = 0
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
                dArea = 0.0
                # Create multipolygon to handle IDL crossing
                if iFlag_verbose:
                    logger.info(f'Feature {i} crosses the international date line, splitting into multiple parts.')
                pMultipolygon = ogr.Geometry(ogr.wkbMultiPolygon)
                aCoord_gcs_split = split_international_date_line_polygon_coordinates(aCoord)
                for aCoord_gcs in aCoord_gcs_split:
                    dArea += calculate_polygon_area(aCoord_gcs[:,0], aCoord_gcs[:,1])
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
                aArea.append(dArea)
                pFeature_mesh = pLayer_mesh.GetNextFeature()
            else:
                dArea = calculate_polygon_area(aCoord[:,0], aCoord[:,1])
                wkt = pPolygon.ExportToWkt()
                aPolygon.append((current_cellid, wkt))
                aArea.append(dArea)
                pFeature_mesh = pLayer_mesh.GetNextFeature()
        else:
            #a feature may still have multipolygon
            dArea = 0.0
            for iPart in range(pPolygon.GetGeometryCount()):
                pPolygon_part = pPolygon.GetGeometryRef(iPart)
                aCoords_part = get_geometry_coordinates(pPolygon_part)
                dArea += calculate_polygon_area(aCoords_part[:,0], aCoords_part[:,1])
            wkt = pPolygon.ExportToWkt()
            aPolygon.append((current_cellid, wkt))
            aArea.append(dArea)
            pFeature_mesh = pLayer_mesh.GetNextFeature()

        i += 1

        # Progress reporting during feature pre-processing
        if i % 1000 == 0:
            logger.info(f"Pre-processed {i} features...")

    if iFlag_verbose:
        logger.info(f"run_remap: Pre-processing completed. Found {len(aPolygon)} features to process")

    return aPolygon, aArea, sProjection_source_wkt

def run_remap(sFilename_target_mesh,
                                sFilename_source_mesh,
                                aFilename_source_raster,
                                dArea_mean,
                                iFlag_remap_method_in = 1,
                            iFlag_stat_in = 1,
                            iFlag_save_clipped_raster_in=0,
                            sFolder_raster_out_in = None,
                            iFlag_verbose=False,
                            iFeature_parallel_threshold=4000):
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


    if os.path.exists(sFilename_target_mesh):
        #remove the file using the vector driver
        pDriver_vector.DeleteDataSource(sFilename_target_mesh)



    #get the raster file extension
    sFilename_raster = aFilename_source_raster[0]  #just use the first raster to get the extension
    sExtension = os.path.splitext(sFilename_raster)[1]
    sName = os.path.basename(sFilename_raster)
    sRasterName_no_extension = os.path.splitext(sName)[0]

    if iFlag_verbose:
        logger.info("run_remap: Reading raster metadata and determining processing bounds...")

    #get the highest resolution raster to determine the pixel size
    dPixelWidth = None
    pPixelHeight = None
    for sFilename_raster in aFilename_source_raster:
        #use sraster class to read the raster info
        pRaster = sraster(sFilename_in=sFilename_raster)
        pRaster.read_metadata()
        if dPixelWidth is None or pRaster.dResolution_x < dPixelWidth:
            dPixelWidth = pRaster.dResolution_x
        if pPixelHeight is None or abs(pRaster.dResolution_y) < abs(pPixelHeight):
            pPixelHeight = pRaster.dResolution_y
        dMissing_value = pRaster.dNoData

    # Determine optimal resampling method based on resolution comparison
    # This will override iFlag_remap_method if mesh resolution is too coarse
    sRemap_method_auto, iRemap_method_auto = _determine_optimal_resampling(dArea_mean, dPixelWidth, abs(pPixelHeight), iFlag_verbose)

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

    aPolygon, aArea, sProjection_source_wkt = get_polygon_list(sFilename_source_mesh, iFlag_verbose)
    pSpatialRef_target = osr.SpatialReference()
    pSpatialRef_target.ImportFromWkt(sProjection_source_wkt)

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


    options = ['COMPRESS=DEFLATE', 'PREDICTOR=2'] #reseverd for future use

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


    logger.info("run_remap: Starting main feature processing loop...")

    #use multiprocessing to speed up the processing
    start_time = time.time()
    successful_features = 0
    failed_features = []

    # Prepare a serializable copy of warp options (convert dstSRS to WKT if needed)
    gdal_warp_options_serial = gdal_warp_options_base.copy()
    if 'dstSRS' in gdal_warp_options_serial and hasattr(gdal_warp_options_serial['dstSRS'], 'ExportToWkt'):
        try:
            gdal_warp_options_serial['dstSRS'] = gdal_warp_options_serial['dstSRS'].ExportToWkt()
        except Exception:
            gdal_warp_options_serial['dstSRS'] = str(gdal_warp_options_serial['dstSRS'])

    n_features = len(aPolygon)
    max_workers = min(cpu_count(), max(1, n_features))
    logger.info(f"Preparing to process {n_features} features (parallel threshold={iFeature_parallel_threshold})")

    # Build ordered task list (keeps original order)
    tasks = []
    for idx, (cellid, wkt) in enumerate(aPolygon):
        tasks.append((idx, cellid, wkt, aFilename_source_raster, gdal_warp_options_serial, dMissing_value, iFlag_verbose))

    # Choose serial or parallel processing based on threshold
    if n_features <= iFeature_parallel_threshold:
        logger.info(f"Feature count ({n_features}) <= threshold ({iFeature_parallel_threshold}); using serial processing")
        for task in tasks:
            feature_idx, cellid, success, payload = _process_task(task)

            if not success:
                failed_features.append({"feature_id": cellid, "error": payload, "envelope": None})
                if iFlag_verbose:
                    logger.warning(f"Feature {cellid} failed: {payload}")
                continue

            # payload is stats dict
            stats = payload
            try:
                # write feature geometry and attributes to output layer
                pFeature_out = ogr.Feature(pLayer_defn_out)
                # set geometry from WKT
                geom = ogr.CreateGeometryFromWkt(aPolygon[feature_idx][1])
                pFeature_out.SetGeometry(geom)
                pFeature_out.SetField('cellid', int(cellid))
                pFeature_out.SetField('area', aArea[feature_idx])
                if iFlag_stat_in == 1:
                    pFeature_out.SetField('mean', float(stats.get('mean', np.nan)))
                    pFeature_out.SetField('min', float(stats.get('min', np.nan)))
                    pFeature_out.SetField('max', float(stats.get('max', np.nan)))
                    pFeature_out.SetField('std', float(stats.get('std', np.nan)))
                pLayer_out.CreateFeature(pFeature_out)
                pFeature_out = None
                successful_features += 1
            except Exception as e:
                failed_features.append({"feature_id": cellid, "error": str(e), "envelope": None})
                logger.error(f"Failed writing feature {cellid}: {e}")
    else:
        logger.info(f"Feature count ({n_features}) > threshold ({iFeature_parallel_threshold}); using multiprocessing with {max_workers} workers")
        # Use ProcessPoolExecutor.map to preserve task order in results
        with ProcessPoolExecutor(max_workers=max_workers) as exe:
            # exe.map will yield results in same order as tasks
            for result in exe.map(_process_task, tasks):
                feature_idx, cellid, success, payload = result

                if not success:
                    failed_features.append({"feature_id": cellid, "error": payload, "envelope": None})
                    if iFlag_verbose:
                        logger.warning(f"Feature {cellid} failed: {payload}")
                    continue

                # payload is stats dict
                stats = payload
                try:
                    # write feature geometry and attributes to output layer
                    pFeature_out = ogr.Feature(pLayer_defn_out)
                    # set geometry from WKT
                    geom = ogr.CreateGeometryFromWkt(aPolygon[feature_idx][1])
                    pFeature_out.SetGeometry(geom)
                    pFeature_out.SetField('cellid', int(cellid))
                    pFeature_out.SetField('area', aArea[feature_idx])
                    if iFlag_stat_in == 1:
                        pFeature_out.SetField('mean', float(stats.get('mean', np.nan)))
                        pFeature_out.SetField('min', float(stats.get('min', np.nan)))
                        pFeature_out.SetField('max', float(stats.get('max', np.nan)))
                        pFeature_out.SetField('std', float(stats.get('std', np.nan)))
                    pLayer_out.CreateFeature(pFeature_out)
                    pFeature_out = None
                    successful_features += 1
                except Exception as e:
                    failed_features.append({"feature_id": cellid, "error": str(e), "envelope": None})
                    logger.error(f"Failed writing feature {cellid}: {e}")

    # end multiprocessing block

    # flush and close output
    pDataset_out.FlushCache()
    pDataset_out = None

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
        base_name = os.path.splitext(sFilename_target_mesh)[0]
        failure_report_file = f"{base_name}_failures.log"
        try:
            with open(failure_report_file, 'w') as f:
                f.write(f"Processing failure report - {time.ctime()}\n")
                f.write(f"Total features processed: {len(aPolygon)}\n")
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

