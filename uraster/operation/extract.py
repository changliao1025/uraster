# Extract module for uraster - contains remap workflow functions
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count
import os
import logging
import time
import traceback
from typing import Optional, Tuple, List, Dict, Any, Union
import numpy as np
from osgeo import gdal, ogr, osr
gdal.UseExceptions()
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename
from uraster.classes.sraster import sraster
from uraster.utility import get_polygon_list, get_unique_values_from_rasters

# Try to import psutil for memory monitoring (optional)
try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False


# Set up logging
logger = logging.getLogger(__name__)
crs = "EPSG:4326"

# Initialize GDAL drivers with error handling
try:
    pDriver_geojson = ogr.GetDriverByName('GeoJSON')
    pDriver_shp = ogr.GetDriverByName('ESRI Shapefile')
    if pDriver_geojson is None or pDriver_shp is None:
        raise RuntimeError("Failed to initialize required GDAL drivers")
except Exception as e:
    logger.error(f"Error initializing GDAL drivers: {e}")
    raise

# Constants for processing thresholds
IDL_LONGITUDE_THRESHOLD = 100  # Degrees - threshold for detecting IDL crossing
WARP_TIMEOUT_SECONDS = 30      # Seconds - timeout for GDAL Warp operations
PROGRESS_REPORT_INTERVAL = 5   # Report progress every N features
MAX_CONSECUTIVE_FAILURES = 10   # Maximum consecutive failures before stopping
HEARTBEAT_INTERVAL = 5          # Seconds between heartbeat logs during long operations




def _determine_optimal_resampling(
    dArea_min: float,
    dPixelWidth: float,
    dPixelHeight: float,
    iFlag_verbose: bool = False,
    dResolution_ratio_threshold: float = 3.0
) -> Tuple[str, int]:
    """
    Determine optimal resampling method based on mesh and raster resolution comparison.

    Compares the characteristic mesh cell size with raster resolution to decide
    whether to use nearest neighbor (when raster is much finer) or weighted
    averaging (when mesh and raster resolutions are comparable).

    Args:
        dArea_mean (float): Mean area of mesh cells in square degrees
        dPixelWidth (float): Raster pixel width in degrees
        dPixelHeight (float): Raster pixel height in degrees (absolute value)
        iFlag_verbose (bool, optional): If True, print detailed progress messages.
            If False, only print error messages. Default is False.
        dResolution_ratio_threshold (float, optional): Threshold ratio for switching
            to weighted averaging. Default is 3.0.

    Returns:
        Tuple[str, int]: (resample_method_string, resample_method_code)
            - resample_method_string: 'near' or 'average'
            - resample_method_code: integer code (1 for nearest, 3 for average)

    Raises:
        ValueError: If input parameters are invalid
        ZeroDivisionError: If raster resolution is zero
    """
    # Input validation
    if not isinstance(dArea_min, (int, float)) or dArea_min <= 0:
        raise ValueError(
            f"dArea_mean must be a positive number, got {dArea_min}")

    if not isinstance(dPixelWidth, (int, float)) or dPixelWidth == 0:
        raise ValueError(
            f"dPixelWidth must be a non-zero number, got {dPixelWidth}")

    if not isinstance(dPixelHeight, (int, float)) or dPixelHeight == 0:
        raise ValueError(
            f"dPixelHeight must be a non-zero number, got {dPixelHeight}")

    if not isinstance(dResolution_ratio_threshold, (int, float)) or dResolution_ratio_threshold <= 0:
        raise ValueError(
            f"dResolution_ratio_threshold must be positive, got {dResolution_ratio_threshold}")

    try:
        # Estimate characteristic mesh cell dimension (approximate square root of mean area)
        dMesh_characteristic_size = np.sqrt(dArea_min)

        # Use the coarser of the two raster dimensions
        dRaster_resolution = max(abs(dPixelWidth), abs(dPixelHeight))

        # Calculate resolution ratio (mesh size / raster size)
        dResolution_ratio = dMesh_characteristic_size / dRaster_resolution

        if iFlag_verbose:
            logger.info("="*60)
            logger.info("Resolution Comparison Analysis:")
            logger.info(
                f"  Raster resolution: {dRaster_resolution:.6f} degrees ({dRaster_resolution*111:.2f} km at equator)")
            logger.info(
                f"  Mean mesh cell size: {dMesh_characteristic_size:.6f} degrees ({dMesh_characteristic_size*111:.2f} km at equator)")
            logger.info(
                f"  Resolution ratio (mesh/raster): {dResolution_ratio:.2f}")
            logger.info(
                f"  Threshold for weighted averaging: {dResolution_ratio_threshold:.2f}")

        # Decision logic
        if dResolution_ratio < dResolution_ratio_threshold:
            # Mesh cells are comparable to or smaller than raster resolution
            # Use weighted averaging to properly capture sub-pixel variations
            recommended_method = 'average'
            recommended_code = 3
            if iFlag_verbose:
                logger.warning(
                    f"Mesh resolution is close to raster resolution (ratio: {dResolution_ratio:.2f})")
                logger.warning(
                    "Switching to WEIGHTED AVERAGING (average) for accuracy")
                logger.warning(
                    "Consider using higher resolution raster data for better results")
        else:
            # Raster is much finer than mesh - nearest neighbor is appropriate
            recommended_method = 'near'
            recommended_code = 1
            if iFlag_verbose:
                logger.info(
                    f"Raster is significantly finer than mesh (ratio: {dResolution_ratio:.2f})")
                logger.info(
                    "Using NEAREST NEIGHBOR resampling (sufficient for this resolution ratio)")

        if iFlag_verbose:
            logger.info("="*60)

        return recommended_method, recommended_code

    except Exception as e:
        logger.error(f"Error in _determine_optimal_resampling: {e}")
        logger.error(f"Traceback: {traceback.format_exc()}")
        # Return safe default
        return 'near', 1


def _process_single_polygon(
    iFeature_idx: int,
    iCellid: Union[int, str],
    sWkt: str,
    aFilename_source_raster: List[str],
    gdal_warp_options_base: Dict[str, Any],
    dMissing_value: float,
    iFlag_discrete_in: bool = False,
    iFlag_verbose_in: bool = False
) -> Tuple[int, Union[int, str], bool, Union[Dict[str, float], str]]:
    """
    Process a single polygon with GDAL Warp operation and calculate statistics.

    Args:
        feature_idx (int): Index of the feature being processed
        cellid (Union[int, str]): Cell identifier for logging and tracking
        wkt (str): Well-Known Text representation of the polygon geometry to clip raster
        aFilename_source_raster (List[str]): List of source raster filenames
        gdal_warp_options_base (Dict[str, Any]): Base GDAL warp options dictionary
        dMissing_value (float): NoData/missing value to exclude from statistics
        iFlag_verbose (bool, optional): If True, print detailed progress messages.
            If False, only print error messages. Default is False.

    Returns:
        Tuple[int, Union[int, str], bool, Union[Dict[str, float], str]]:
            (feature_idx, cellid, success_flag, stats_dict_or_error_message)
            - feature_idx: Original feature index
            - cellid: Original cell identifier
            - success_flag: True if processing succeeded, False otherwise
            - stats_dict_or_error_message: Statistics dictionary on success, error message on failure

    Raises:
        ValueError: If input parameters are invalid
    """
    # Input validation
    if not isinstance(sWkt, str) or not sWkt.strip():
        return iFeature_idx, iCellid, False, "Invalid WKT string provided"

    if not isinstance(aFilename_source_raster, list) or not aFilename_source_raster:
        return iFeature_idx, iCellid, False, "Invalid or empty raster filename list"

    if not isinstance(gdal_warp_options_base, dict):
        return iFeature_idx, iCellid, False, "Invalid GDAL warp options dictionary"

    # Initialize cleanup variables
    pPolygonWKT_file = None
    srs_wgs84 = None
    pFeature_clip = None
    pLayer_clip = None
    pDataset_clip = None
    pDataset_warp = None
    polygon = None

    try:
        # Create spatial reference system
        srs_wgs84 = osr.SpatialReference()
        srs_wgs84.ImportFromEPSG(4326)

        # Create geometry from WKT with validation
        polygon = ogr.CreateGeometryFromWkt(sWkt)
        if polygon is None:
            return iFeature_idx, iCellid, False, f"Failed to create geometry from WKT for feature {iCellid}"

        if not polygon.IsValid():
            logger.warning(
                f"Invalid geometry for feature {iCellid}, attempting to fix")
            polygon = polygon.Buffer(0)  # Attempt to fix invalid geometry
            if not polygon.IsValid():
                return iFeature_idx, iCellid, False, f"Cannot fix invalid geometry for feature {iCellid}"

        # Make a copy of the warp options to modify
        gdal_warp_options = gdal_warp_options_base.copy()

        # Create temporary shapefile in memory
        pPolygonWKT_file = f'/vsimem/polygon_wkt_{iCellid}_{iFeature_idx}.shp'
        pDataset_clip = pDriver_shp.CreateDataSource(pPolygonWKT_file)
        if pDataset_clip is None:
            return iFeature_idx, iCellid, False, f"Failed to create temporary dataset for feature {iCellid}"

        pLayer_clip = pDataset_clip.CreateLayer(
            'polygon', geom_type=ogr.wkbPolygon, srs=srs_wgs84)
        if pLayer_clip is None:
            return iFeature_idx, iCellid, False, f"Failed to create layer for feature {iCellid}"

        pFeature_clip = ogr.Feature(pLayer_clip.GetLayerDefn())
        pFeature_clip.SetGeometry(polygon)
        pLayer_clip.CreateFeature(pFeature_clip)
        pDataset_clip.FlushCache()

        gdal_warp_options['cutlineDSName'] = pPolygonWKT_file

        # Run GDAL Warp with timing
        warp_start_time = time.time()
        pWrapOption = gdal.WarpOptions(**gdal_warp_options)
        pDataset_warp = gdal.Warp(
            "", aFilename_source_raster, options=pWrapOption)

        if pDataset_warp is None:
            return iFeature_idx, iCellid, False, f"GDAL Warp failed for feature {iCellid}"

        aData_clip = pDataset_warp.ReadAsArray()
        warp_duration = time.time() - warp_start_time

        if iFlag_verbose_in:
            logger.info(
                f"GDAL Warp completed for feature {iCellid} in {warp_duration:.2f} seconds")

        if aData_clip is None:
            return iFeature_idx, iCellid, False, f"GDAL Warp returned no data for feature {iCellid}"

        # Check for reasonable data dimensions
        if aData_clip.size == 0:
            logger.warning(f"Empty data array for feature {iCellid}")
            return iFeature_idx, iCellid, False, f"Empty data array for feature {iCellid}"

        # Make a copy of the data to allow immediate dataset cleanup
        aData_clip_copy = aData_clip.copy()

        # Calculate statistics with improved error handling
        try:
            # Filter out missing/nodata values with proper handling of different data types
            if np.isnan(dMissing_value):
                valid_mask = ~np.isnan(aData_clip_copy)
            else:
                valid_mask = aData_clip_copy != dMissing_value

            valid_data = aData_clip_copy[valid_mask]

            if len(valid_data) == 0:
                if iFlag_discrete_in == 1:
                    # Handle case for discrete data with no valid pixels
                    stats = {
                        'mode': float(np.nan),
                        'count': 0
                    }
                else:
                    # No valid pixels for this feature: not treated as an error
                    stats = {
                        'mean': float(np.nan),
                        'min': float(np.nan),
                        'max': float(np.nan),
                        'std': float(np.nan),
                        'count': 0
                    }
                if iFlag_verbose_in:
                    logger.info(f"No valid data found for feature {iCellid}")
            else:
                # Compute statistics on valid data with error handling
                if iFlag_discrete_in :
                    # Combine mode calculation and percentage calculation for discrete data
                    values, counts = np.unique(valid_data, return_counts=True)
                    mode_index = np.argmax(counts)
                    stats = {
                        'mode': values[mode_index],
                        'count': int(len(valid_data)),
                        **{f'percentage_{val}': float(cnt) / len(valid_data) * 100.0 for val, cnt in zip(values, counts)}
                    }
                    #also calculate the percentage of each unique value
                    for val, cnt in zip(values, counts):
                        stats[f'percentage_{val}'] = float(cnt) / len(valid_data) * 100.0
                else:
                    stats = {
                    'mean': float(np.mean(valid_data)),
                    'min': float(np.min(valid_data)),
                    'max': float(np.max(valid_data)),
                    'std': float(np.std(valid_data)),
                    'count': int(len(valid_data))
                    }

                if iFlag_verbose_in:
                    logger.info(
                        f"Computed stats for feature {iCellid}: {stats['count']} valid pixels")

        except Exception as stats_error:
            logger.error(
                f"Error computing statistics for feature {iCellid}: {stats_error}")
            return iFeature_idx, iCellid, False, f"Statistics computation failed: {str(stats_error)}"

        return iFeature_idx, iCellid, True, stats

    except Exception as e:
        error_msg = f"Error processing feature {iCellid}: {str(e)}"
        logger.error(error_msg)
        logger.error(f"Traceback: {traceback.format_exc()}")
        return iFeature_idx, iCellid, False, error_msg

    finally:
        # Comprehensive cleanup with error handling
        try:
            if pDataset_warp is not None:
                pDataset_warp = None
        except Exception as e:
            logger.warning(
                f"Error cleaning up warp dataset for feature {iCellid}: {e}")

        try:
            if pFeature_clip is not None:
                pFeature_clip = None
        except Exception as e:
            logger.warning(
                f"Error cleaning up feature for feature {iCellid}: {e}")

        try:
            if pLayer_clip is not None:
                pLayer_clip = None
        except Exception as e:
            logger.warning(
                f"Error cleaning up layer for feature {iCellid}: {e}")

        try:
            if pDataset_clip is not None:
                pDataset_clip = None
        except Exception as e:
            logger.warning(
                f"Error cleaning up dataset for feature {iCellid}: {e}")

        try:
            if polygon is not None:
                polygon = None
        except Exception as e:
            logger.warning(
                f"Error cleaning up geometry for feature {iCellid}: {e}")

        try:
            if srs_wgs84 is not None:
                srs_wgs84 = None
        except Exception as e:
            logger.warning(
                f"Error cleaning up spatial reference for feature {iCellid}: {e}")

        # Clean up temporary file
        try:
            if pPolygonWKT_file is not None:
                gdal.Unlink(pPolygonWKT_file)
        except Exception as e:
            logger.warning(
                f"Error cleaning up temporary file {pPolygonWKT_file}: {e}")


def _merge_raster_parts(
    data_arrays: List[np.ndarray],
    transforms: List[Any],
    feature_id: Union[int, str],
    dMissing_value: float,
    iFlag_verbose: bool = False
) -> Tuple[Optional[np.ndarray], Optional[Any]]:
    """
    Merge multiple raster arrays from IDL-split polygons into a single data array.

    Extracts all valid (non-missing) data values from multiple raster arrays
    and concatenates them into a single 1D array for statistics calculation.

    Args:
        data_arrays (List[np.ndarray]): List of numpy arrays from different polygon parts
        transforms (List[Any]): List of geotransforms (one per data array)
        feature_id (Union[int, str]): Feature identifier for logging
        dMissing_value (float): Missing/NoData value to exclude from results
        iFlag_verbose (bool, optional): If True, print detailed progress messages.
            If False, only print error messages. Default is False.

    Returns:
        Tuple[Optional[np.ndarray], Optional[Any]]: (merged_1D_array, first_transform)
        - merged_1D_array: 1D array containing all valid data values, None on failure
        - first_transform: Geotransform from first polygon part (dummy value), None on failure

    Raises:
        ValueError: If input parameters are invalid
    """
    # Input validation
    if not isinstance(data_arrays, list) or not data_arrays:
        logger.error(f"Invalid or empty data_arrays for feature {feature_id}")
        return None, None

    if not isinstance(transforms, list) or len(transforms) != len(data_arrays):
        logger.error(
            f"Transforms list length mismatch for feature {feature_id}")
        return None, None

    try:
        # Handle single array case
        if len(data_arrays) == 1:
            data_array = data_arrays[0]

            if not isinstance(data_array, np.ndarray):
                logger.error(
                    f"Invalid data array type for feature {feature_id}")
                return None, None

            if data_array.size == 0:
                logger.warning(f"Empty data array for feature {feature_id}")
                return np.array([]), transforms[0] if transforms else None

            # Extract valid data with proper handling of NaN missing values
            if np.isnan(dMissing_value):
                valid_mask = ~np.isnan(data_array)
            else:
                valid_mask = data_array != dMissing_value

            valid_data = data_array[valid_mask]

            if valid_data.size > 0:
                result_data = valid_data.flatten()
            else:
                result_data = np.array([])

            if iFlag_verbose:
                logger.info(
                    f"Single array processed for feature {feature_id}: {result_data.size} valid pixels")

            return result_data, transforms[0] if transforms else None

        # Handle multiple arrays case
        all_valid_data = []
        total_pixels = 0

        for i, data_array in enumerate(data_arrays):
            if not isinstance(data_array, np.ndarray):
                logger.warning(
                    f"Skipping invalid data array {i} for feature {feature_id}")
                continue

            if data_array.size == 0:
                logger.warning(
                    f"Skipping empty data array {i} for feature {feature_id}")
                continue

            # Extract all valid (non-nodata) values with proper NaN handling
            try:
                if np.isnan(dMissing_value):
                    valid_mask = ~np.isnan(data_array)
                else:
                    valid_mask = data_array != dMissing_value

                valid_data = data_array[valid_mask]

                if valid_data.size > 0:
                    all_valid_data.append(valid_data.flatten())
                    total_pixels += valid_data.size

                    if iFlag_verbose:
                        logger.debug(
                            f"Array {i} for feature {feature_id}: {valid_data.size} valid pixels")

            except Exception as array_error:
                logger.warning(
                    f"Error processing array {i} for feature {feature_id}: {array_error}")
                continue

        if not all_valid_data:
            logger.warning(
                f"No valid data found in any part for IDL feature {feature_id}")
            return np.array([]), transforms[0] if transforms else None

        # Concatenate all valid data into a single 1D array
        try:
            merged_data = np.concatenate(all_valid_data)

            if iFlag_verbose:
                logger.info(f"Successfully merged {len(data_arrays)} raster parts for feature {feature_id}: "
                            f"{merged_data.size} valid pixels from {total_pixels} total pixels")

            return merged_data, transforms[0] if transforms else None

        except Exception as concat_error:
            logger.error(
                f"Error concatenating arrays for feature {feature_id}: {concat_error}")
            return None, None

    except Exception as e:
        logger.error(
            f"Error merging raster parts for feature {feature_id}: {str(e)}")
        logger.error(f"Traceback: {traceback.format_exc()}")
        return None, None


def _process_multipolygon_idl(
    iFeature_idx: int,
    iCellid: Union[int, str],
    sWkt: str,
    aFilename_source_raster: List[str],
    gdal_warp_options_base: Dict[str, Any],
    dMissing_value: float,
    iFlag_discrete_in: bool = False,
    iFlag_verbose_in: bool = False
) -> Tuple[Optional[np.ndarray], Optional[Any]]:
    """
    Process a multipolygon that crosses the International Date Line (IDL).

    Handles IDL-crossing features by processing each polygon part separately
    and merging the results into a single data array for statistics calculation.

    Args:
        feature_idx (int): Index of the feature being processed
        feature_id (Union[int, str]): Feature identifier for logging
        wkt (str): Well-Known Text representation of the multipolygon geometry
        aFilename_source_raster (List[str]): List of source raster filenames
        gdal_warp_options_base (Dict[str, Any]): Base GDAL warp options dictionary
        dMissing_value (float): Missing/NoData value for the raster
        iFlag_verbose (bool, optional): If True, print detailed progress messages.
            If False, only print error messages. Default is False.

    Returns:
        Tuple[Optional[np.ndarray], Optional[Any]]: (merged_data_array, merged_geotransform)
        - merged_data_array: 1D array containing all valid data values, None on failure
        - merged_geotransform: Geotransform from first polygon part (dummy value), None on failure

    Raises:
        ValueError: If input parameters are invalid
    """
    # Input validation
    if not isinstance(sWkt, str) or not sWkt.strip():
        logger.error(
            f"Invalid WKT string for multipolygon feature {iCellid}")
        return None, None

    if not isinstance(aFilename_source_raster, list) or not aFilename_source_raster:
        logger.error(
            f"Invalid raster filename list for multipolygon feature {iCellid}")
        return None, None

    multipolygon = None
    try:
        # Create and validate multipolygon geometry
        multipolygon = ogr.CreateGeometryFromWkt(sWkt)
        if multipolygon is None:
            logger.error(
                f"Failed to create multipolygon geometry from WKT for feature {iCellid}")
            return None, None

        if not multipolygon.IsValid():
            logger.warning(
                f"Invalid multipolygon geometry for feature {iCellid}, attempting to fix")
            multipolygon = multipolygon.Buffer(0)  # Attempt to fix
            if not multipolygon.IsValid():
                logger.error(
                    f"Cannot fix invalid multipolygon geometry for feature {iCellid}")
                return None, None

        nGeometries = multipolygon.GetGeometryCount()
        if nGeometries == 0:
            logger.warning(
                f"Multipolygon has no geometry parts for feature {iCellid}")
            return None, None

        if iFlag_verbose_in:
            logger.info(
                f"Processing {nGeometries} polygon parts for IDL-crossing feature {iCellid}")

        merged_data_arrays = []
        merged_transforms = []
        successful_parts = 0

        # Process each polygon part separately
        for iPart in range(nGeometries):
            try:
                polygon_part = multipolygon.GetGeometryRef(iPart)
                if polygon_part is None:
                    logger.warning(
                        f"Polygon part {iPart} is None for feature {iCellid}")
                    continue

                if not polygon_part.IsValid():
                    logger.warning(
                        f"Invalid polygon part {iPart} for feature {iCellid}, attempting to fix")
                    polygon_part = polygon_part.Buffer(0)
                    if not polygon_part.IsValid():
                        logger.warning(
                            f"Cannot fix polygon part {iPart} for feature {iCellid}, skipping")
                        continue

                # Process this polygon part - convert geometry to WKT
                polygon_part_wkt = polygon_part.ExportToWkt()
                if not polygon_part_wkt:
                    logger.warning(
                        f"Failed to export WKT for polygon part {iPart} of feature {iCellid}")
                    continue

                part_result = _process_single_polygon(
                    iFeature_idx, f"{iCellid}_part{iPart}", polygon_part_wkt, aFilename_source_raster,
                      gdal_warp_options_base,
                    dMissing_value, iFlag_discrete_in=iFlag_discrete_in, iFlag_verbose_in=iFlag_verbose_in)


                if len(part_result) != 4 or not part_result[2]:
                    logger.warning(
                        f"Failed to process polygon part {iPart} of feature {iCellid}: {part_result[3] if len(part_result) > 3 else 'Unknown error'}")
                    continue

                # Extract the stats from the result and convert to data array for merging
                part_stats = part_result[3]
                if isinstance(part_stats, dict) and 'count' in part_stats and 'mean' in part_stats:
                    if part_stats['count'] > 0 and not np.isnan(part_stats['mean']):
                        # Create a dummy array with the mean value repeated 'count' times
                        # This preserves the statistical properties for merging
                        part_data = np.full(
                            part_stats['count'], part_stats['mean'])
                        merged_data_arrays.append(part_data)
                        # Transform not needed for statistics
                        merged_transforms.append(None)
                        successful_parts += 1

                        if iFlag_verbose_in:
                            logger.debug(
                                f"Successfully processed part {iPart} of feature {iCellid}: {part_stats['count']} pixels")
                    else:
                        if iFlag_verbose_in:
                            logger.debug(
                                f"Part {iPart} of feature {iCellid} has no valid data")
                else:
                    logger.warning(
                        f"Invalid stats returned for polygon part {iPart} of feature {iCellid}: {type(part_stats)}")
                    continue

            except Exception as part_error:
                logger.warning(
                    f"Error processing polygon part {iPart} of feature {iCellid}: {part_error}")
                continue

        if not merged_data_arrays:
            logger.warning(
                f"No polygon parts could be processed for IDL feature {iCellid}")
            return np.array([]), None

        if iFlag_verbose:
            logger.info(
                f"Successfully processed {successful_parts}/{nGeometries} parts for IDL feature {iCellid}")

        # Merge the data arrays and transforms
        merged_data, merged_transform = _merge_raster_parts(
            merged_data_arrays, merged_transforms, iCellid, dMissing_value, iFlag_verbose)

        return merged_data, merged_transform

    except Exception as e:
        logger.error(
            f"Error processing multipolygon IDL feature {iCellid}: {str(e)}")
        logger.error(f"Traceback: {traceback.format_exc()}")
        return None, None

    finally:
        # Clean up geometry
        try:
            if multipolygon is not None:
                multipolygon = None
        except Exception as e:
            logger.warning(
                f"Error cleaning up multipolygon geometry for feature {iCellid}: {e}")


def _process_task(args: Tuple[int, Union[int, str], str, List[str], Dict[str, Any], float, bool, bool]) -> Tuple[int, Union[int, str], bool, Union[Dict[str, float], str]]:
    """
    Module-level worker function for multiprocessing polygon processing.

    Processes a single polygon or multipolygon feature by determining its geometry type
    and calling the appropriate processing function. This function is designed to be
    used with multiprocessing pools.

    Args:
        args (Tuple): Tuple containing:
            - feature_idx (int): Index of the feature being processed
            - cellid (Union[int, str]): Cell identifier for the feature
            - wkt (str): Well-Known Text representation of the geometry
            - aFilename_source_raster (List[str]): List of source raster files
            - gdal_warp_options_base (Dict[str, Any]): GDAL warp options dictionary
            - dMissing_value (float): Missing/NoData value for the raster
            - iFlag_verbose (bool): Verbose logging flag

    Returns:
        Tuple[int, Union[int, str], bool, Union[Dict[str, float], str]]:
            (feature_idx, cellid, success_flag, stats_dict_or_error_message)
            - feature_idx: Original feature index
            - cellid: Original cell identifier
            - success_flag: True if processing succeeded, False otherwise
            - stats_dict_or_error_message: Statistics dictionary on success, error message on failure

    Note:
        This function handles both POLYGON and MULTIPOLYGON geometries, with special
        handling for International Date Line crossing multipolygons.
    """
    # Unpack arguments with validation
    try:
        feature_idx, cellid, wkt, aFilename_source_raster, gdal_warp_options_base, dMissing_value, iFlag_discrete_in, iFlag_verbose = args
    except (ValueError, TypeError) as e:
        logger.error(f"Invalid arguments passed to _process_task: {e}")
        return -1, "unknown", False, f"Invalid arguments: {str(e)}"

    # Input validation
    if not isinstance(wkt, str) or not wkt.strip():
        return feature_idx, cellid, False, "Invalid WKT string provided"

    if not isinstance(aFilename_source_raster, list) or not aFilename_source_raster:
        return feature_idx, cellid, False, "Invalid or empty raster filename list"

    pPolygon = None
    try:
        # Create and validate geometry
        pPolygon = ogr.CreateGeometryFromWkt(wkt)
        if pPolygon is None:
            return feature_idx, cellid, False, f"Failed to create geometry from WKT for feature {cellid}"

        if not pPolygon.IsValid():
            logger.warning(
                f"Invalid geometry for feature {cellid}, attempting to fix")
            pPolygon = pPolygon.Buffer(0)  # Attempt to fix invalid geometry
            if not pPolygon.IsValid():
                return feature_idx, cellid, False, f"Cannot fix invalid geometry for feature {cellid}"

        sGeometry_type = pPolygon.GetGeometryName()

        if iFlag_verbose:
            logger.debug(
                f"Processing {sGeometry_type} geometry for feature {cellid}")

        if sGeometry_type == "POLYGON":
            return _process_single_polygon(feature_idx, cellid,
                                           wkt, aFilename_source_raster, gdal_warp_options_base, dMissing_value, iFlag_discrete_in, iFlag_verbose)

        elif sGeometry_type == "MULTIPOLYGON":
            merged_data, merged_transform = _process_multipolygon_idl(feature_idx, cellid,
                                                                      wkt, aFilename_source_raster, gdal_warp_options_base, dMissing_value, iFlag_discrete_in, iFlag_verbose)

            if merged_data is None:
                return feature_idx, cellid, False, "Failed to process multipolygon"

            # Calculate statistics for multipolygon data with improved error handling
            try:
                # Handle different missing value types
                if np.isnan(dMissing_value):
                    valid_mask = ~np.isnan(merged_data)
                else:
                    valid_mask = merged_data != dMissing_value

                valid_data = merged_data[valid_mask]

                if len(valid_data) == 0:
                    return feature_idx, cellid, False, "No valid data found in multipolygon"

                stats = {
                    'mean': float(np.mean(valid_data)),
                    'min': float(np.min(valid_data)),
                    'max': float(np.max(valid_data)),
                    'std': float(np.std(valid_data)),
                    'count': int(len(valid_data))
                }

                if iFlag_verbose:
                    logger.debug(
                        f"Computed multipolygon stats for feature {cellid}: {stats['count']} valid pixels")

                return feature_idx, cellid, True, stats

            except Exception as stats_error:
                logger.error(
                    f"Error computing multipolygon statistics for feature {cellid}: {stats_error}")
                return feature_idx, cellid, False, f"Statistics computation failed: {str(stats_error)}"
        else:
            logger.error(
                f"Unsupported geometry type for feature {cellid}: {sGeometry_type}")
            return feature_idx, cellid, False, f"Unsupported geometry type: {sGeometry_type}"

    except TimeoutError as e:
        logger.error(f"Timeout processing feature {cellid}: {e}")
        return feature_idx, cellid, False, f"Timeout: {str(e)}"
    except Exception as e:
        logger.error(f"Error processing feature {cellid}: {e}")
        logger.error(f"Traceback: {traceback.format_exc()}")
        return feature_idx, cellid, False, f"Error: {str(e)}"

    finally:
        # Clean up geometry
        try:
            if pPolygon is not None:
                pPolygon = None
        except Exception as e:
            logger.warning(
                f"Error cleaning up geometry for feature {cellid}: {e}")



def run_remap(sFilename_target_mesh,
              sFilename_source_mesh,
              aFilename_source_raster,
              dArea_min,
              iFlag_remap_method_in = 1,
              iFlag_stat_in = True,
              iFlag_save_clipped_raster_in = 0,
              sFolder_raster_out_in = None,
              iFlag_discrete_in = False,
              iFlag_verbose_in = False,
              iFeature_parallel_threshold = 5000):
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
        iFlag_stat_in (bool, optional): Flag to compute statistics (True=yes, False=no).
            Default is True.
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
        logger.error(
            "Invalid remap method specified. Must be 1 (nearest), 2 (bilinear), or 3 (weighted average).")
        return
    iFlag_remap_method = iFlag_remap_method_in

    if iFlag_verbose_in  :
        logger.info("run_remap: Starting input file validation...")
    # check input files
    for idx, sFilename_raster in enumerate(aFilename_source_raster):
        if iFlag_verbose_in  :
            logger.info(
                f"Checking raster file {idx+1}/{len(aFilename_source_raster)}: {os.path.basename(sFilename_raster)}")
        if os.path.exists(sFilename_raster):
            pass
        else:
            logger.error('The raster file does not exist!')
            return

    if iFlag_verbose_in  :
        logger.info(
            f"Checking source mesh file: {os.path.basename(sFilename_source_mesh)}")
    if os.path.exists(sFilename_source_mesh):
        pass
    else:
        logger.error('The vector mesh file does not exist!')
        return
    if iFlag_verbose_in  :
        logger.info("Input file validation completed successfully")

    # Determine output vector format from filename extension
    pDriver_vector = get_vector_driver_from_filename(sFilename_target_mesh)

    if os.path.exists(sFilename_target_mesh):
        # remove the file using the vector driver
        pDriver_vector.DeleteDataSource(sFilename_target_mesh)

    # get the raster file extension
    # just use the first raster to get the extension
    sFilename_raster = aFilename_source_raster[0]
    sExtension = os.path.splitext(sFilename_raster)[1]
    sName = os.path.basename(sFilename_raster)
    sRasterName_no_extension = os.path.splitext(sName)[0]

    if iFlag_verbose_in  :
        logger.info(
            "run_remap: Reading raster metadata and determining processing bounds...")

    # get the highest resolution raster to determine the pixel size
    dPixelWidth = None
    pPixelHeight = None
    for sFilename_raster in aFilename_source_raster:
        # use sraster class to read the raster info
        pRaster = sraster(sFilename_in=sFilename_raster)
        pRaster.read_metadata()
        if dPixelWidth is None or pRaster.dResolution_x < dPixelWidth:
            dPixelWidth = pRaster.dResolution_x
        if pPixelHeight is None or abs(pRaster.dResolution_y) < abs(pPixelHeight):
            pPixelHeight = pRaster.dResolution_y
        dMissing_value = pRaster.dNoData

    # Determine optimal resampling method based on resolution comparison
    # This will override iFlag_remap_method if mesh resolution is too coarse
    sRemap_method_auto, iRemap_method_auto = _determine_optimal_resampling(
        dArea_min, dPixelWidth, abs(pPixelHeight), iFlag_verbose_in  )

    # Use automatically determined method if it's more conservative than user setting
    # Priority: weighted averaging > nearest neighbor
    if iRemap_method_auto == 3 and iFlag_remap_method != 3:
        logger.warning(
            f"Overriding user remap method ({iFlag_remap_method}) with automatic selection (3 - weighted average)")
        logger.warning(
            "This is necessary due to mesh/raster resolution compatibility")
        sRemap_method = sRemap_method_auto
    else:
        # Use user's preferred method
        if iFlag_remap_method == 1:
            sRemap_method = 'near'
        elif iFlag_remap_method == 2:
            sRemap_method = 'near'
        elif iFlag_remap_method == 3:
            sRemap_method = 'average'
        if iFlag_verbose_in  :
            logger.info(f"Using user-specified remap method: {sRemap_method}")

    if iFlag_verbose_in  :
        logger.info("run_remap: Opening mesh dataset and analyzing features...")

    aPolygon, aArea, sProjection_source_wkt = get_polygon_list(
        sFilename_source_mesh, iFlag_verbose_in  )
    pSpatialRef_target = osr.SpatialReference()
    pSpatialRef_target.ImportFromWkt(sProjection_source_wkt)

    # create a polygon feature to save the output
    pDataset_out = pDriver_vector.CreateDataSource(sFilename_target_mesh)
    pLayer_out = pDataset_out.CreateLayer(
        'uraster', pSpatialRef_target, ogr.wkbPolygon)
    pLayer_defn_out = pLayer_out.GetLayerDefn()
    pFeature_out = ogr.Feature(pLayer_defn_out)

    # add id, area and mean, min, max, std of the raster
    pLayer_out.CreateField(ogr.FieldDefn('cellid', ogr.OFTInteger))
    # define a field
    pField = ogr.FieldDefn('area', ogr.OFTReal)
    pField.SetWidth(32)
    pField.SetPrecision(2)
    pLayer_out.CreateField(pField)

    # in the future, we will also copy other attributes from the input geojson file
    if iFlag_discrete_in:
        #we might need to get the unique values first to create the fields
        unique_values = get_unique_values_from_rasters(aFilename_source_raster, dMissing_value, iFlag_verbose_in )
        nValues = len(unique_values)
        logger.info(f"Found {nValues} unique values in raster")
        pLayer_out.CreateField(ogr.FieldDefn('class', ogr.OFTInteger))
        #how about also adding percentage fields for each unique value?
        for val in unique_values:
            field_name = f'frac_{int(val)}'
            pLayer_out.CreateField(ogr.FieldDefn(field_name, ogr.OFTReal))
    else:
        if iFlag_stat_in:
            pLayer_out.CreateField(ogr.FieldDefn('mean', ogr.OFTReal))
            pLayer_out.CreateField(ogr.FieldDefn('min', ogr.OFTReal))
            pLayer_out.CreateField(ogr.FieldDefn('max', ogr.OFTReal))
            pLayer_out.CreateField(ogr.FieldDefn('std', ogr.OFTReal))
        else:
            pLayer_out.CreateField(ogr.FieldDefn('mean', ogr.OFTReal))
            pass

    options = ['COMPRESS=DEFLATE', 'PREDICTOR=2']  # reseverd for future use

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

    # use multiprocessing to speed up the processing
    start_time = time.time()
    successful_features = 0
    failed_features = []

    # Prepare a serializable copy of warp options (convert dstSRS to WKT if needed)
    gdal_warp_options_serial = gdal_warp_options_base.copy()
    if 'dstSRS' in gdal_warp_options_serial and hasattr(gdal_warp_options_serial['dstSRS'], 'ExportToWkt'):
        try:
            gdal_warp_options_serial['dstSRS'] = gdal_warp_options_serial['dstSRS'].ExportToWkt(
            )
        except Exception:
            gdal_warp_options_serial['dstSRS'] = str(
                gdal_warp_options_serial['dstSRS'])

    n_features = len(aPolygon)
    max_workers = min(cpu_count(), max(1, n_features))
    logger.info(
        f"Preparing to process {n_features} features (parallel threshold={iFeature_parallel_threshold})")

    # Build ordered task list (keeps original order)
    tasks = []
    for idx, (cellid, wkt) in enumerate(aPolygon):
        tasks.append((idx, cellid, wkt, aFilename_source_raster, gdal_warp_options_serial, dMissing_value, iFlag_discrete_in, iFlag_verbose_in))

    # Choose serial or parallel processing based on threshold
    if n_features <= iFeature_parallel_threshold:
        logger.info(
            f"Feature count ({n_features}) <= threshold ({iFeature_parallel_threshold}); using serial processing")
        for task in tasks:
            feature_idx, cellid, success, payload = _process_task(task)
            if not success:
                failed_features.append(
                    {"feature_id": cellid, "error": payload, "envelope": None})
                if iFlag_verbose_in:
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
                if iFlag_discrete_in:
                    # populate class and fraction fields
                    pass
                else:
                    if iFlag_stat_in:
                        pFeature_out.SetField(
                            'mean', float(stats.get('mean', np.nan)))
                        pFeature_out.SetField(
                            'min', float(stats.get('min', np.nan)))
                        pFeature_out.SetField(
                            'max', float(stats.get('max', np.nan)))
                        pFeature_out.SetField(
                            'std', float(stats.get('std', np.nan)))
                pLayer_out.CreateFeature(pFeature_out)
                pFeature_out = None
                successful_features += 1
            except Exception as e:
                failed_features.append(
                    {"feature_id": cellid, "error": str(e), "envelope": None})
                logger.error(f"Failed writing feature {cellid}: {e}")
    else:
        logger.info(
            f"Feature count ({n_features}) > threshold ({iFeature_parallel_threshold}); using multiprocessing with {max_workers} workers")
        # Use ProcessPoolExecutor.map to preserve task order in results
        with ProcessPoolExecutor(max_workers=max_workers) as exe:
            # exe.map will yield results in same order as tasks
            for result in exe.map(_process_task, tasks):
                feature_idx, cellid, success, payload = result

                if not success:
                    failed_features.append(
                        {"feature_id": cellid, "error": payload, "envelope": None})
                    if iFlag_verbose_in  :
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
                    if iFlag_stat_in:
                        pFeature_out.SetField(
                            'mean', float(stats.get('mean', np.nan)))
                        pFeature_out.SetField(
                            'min', float(stats.get('min', np.nan)))
                        pFeature_out.SetField(
                            'max', float(stats.get('max', np.nan)))
                        pFeature_out.SetField(
                            'std', float(stats.get('std', np.nan)))
                    pLayer_out.CreateFeature(pFeature_out)
                    pFeature_out = None
                    successful_features += 1
                except Exception as e:
                    failed_features.append(
                        {"feature_id": cellid, "error": str(e), "envelope": None})
                    logger.error(f"Failed writing feature {cellid}: {e}")

    # end multiprocessing block

    # flush and close output
    pDataset_out.FlushCache()
    pDataset_out = None

    # Clean up spatial reference objects to prevent memory leaks
    pSpatialRef_target = None

    # Report processing summary
    total_time = time.time() - start_time
    if iFlag_verbose_in  :
        logger.info(f"Processing completed in {total_time:.2f} seconds")
        logger.info(f"Successfully processed: {successful_features} features")
        logger.info(f"Failed features: {len(failed_features)}")

    if failed_features:
        if iFlag_verbose_in  :
            logger.warning("Failed features summary:")
            for failed in failed_features[:10]:  # Show first 10 failures
                logger.warning(
                    f"  Feature {failed['feature_id']}: {failed['error']}")
            if len(failed_features) > 10:
                logger.warning(
                    f"  ... and {len(failed_features) - 10} more failures")

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
                    f.write(
                        f"Feature {failed['feature_id']}: {failed['error']}\n")
                    f.write(f"  Envelope: {failed['envelope']}\n\n")
            if iFlag_verbose_in  :
                logger.info(f"Failure report saved to: {failure_report_file}")
        except Exception as e:
            logger.error(f"Could not save failure report: {e}")

    return
