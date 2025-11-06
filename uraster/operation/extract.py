
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

def _determine_optimal_resampling(dArea_mean, dPixelWidth, dPixelHeight, iFlag_verbose=False):
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

def _process_task(args):
    """
    Module-level worker for multiprocessing. Accepts a tuple:
    (feature_idx, cellid, wkt, aFilename_source_raster, gdal_warp_options_serial, dMissing_value, iFlag_verbose)
    Returns a tuple (feature_idx, cellid, success, stats_dict_or_error)
    """
    try:
        (feature_idx, cellid, wkt, aFilename_source_raster,
         gdal_warp_options_serial, dMissing_value, iFlag_verbose) = args

        # Recreate geometry from WKT in the worker process
        polygon = ogr.CreateGeometryFromWkt(wkt)
        if polygon is None:
            return (feature_idx, cellid, False, "Invalid geometry WKT")

        # Call existing processing routine (returns array, geotransform) or (None, None) on failure
        aData_clip, newGeoTransform = _process_single_polygon(
            polygon, aFilename_source_raster, gdal_warp_options_serial, feature_idx, iFlag_verbose=iFlag_verbose
        )

        if aData_clip is None:
            return (feature_idx, cellid, False, "Warp failed or returned empty array")

        # Convert to numpy array and mask missing values if provided
        arr = np.asarray(aData_clip)
        if dMissing_value is not None:
            # replace nodata value with np.nan for robust stats
            try:
                arr = np.where(arr == dMissing_value, np.nan, arr)
            except Exception:
                pass

        # flatten and compute stats using nan-aware functions
        flat = arr.flatten()
        valid = flat[~np.isnan(flat)]
        if valid.size == 0:
            stats = {"mean": np.nan, "min": np.nan, "max": np.nan, "std": np.nan, "count": 0}
        else:
            stats = {
                "mean": float(np.nanmean(flat)),
                "min": float(np.nanmin(flat)),
                "max": float(np.nanmax(flat)),
                "std": float(np.nanstd(flat, ddof=0)),
                "count": int(valid.size),
            }

        return (feature_idx, cellid, True, stats)

    except Exception as e:
        return (feature_idx if 'feature_idx' in locals() else -1,
                cellid if 'cellid' in locals() else None,
                False, str(e))

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
                dArea_mean,
                iFlag_remap_method_in = 1,
              iFlag_stat_in = 1,
              iFlag_save_clipped_raster_in=0,
              sFolder_raster_out_in = None,
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

    #get the 


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


    aPolygon = get_polygon_list(sFilename_source_mesh, iFlag_verbose)

    options = ['COMPRESS=DEFLATE', 'PREDICTOR=2']


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

    if iFlag_verbose:
        logger.info("run_remap: Pre-fetching features and analyzing geometries...")

    start_time = time.time()

    # Add crash tracking variables
    failed_features = []
    successful_features = 0
    logger.info("run_remap: Starting main feature processing loop...")

    #use multiprocessing to speed up the processing






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

