from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import Pool, cpu_count
import os
import logging
import time
import traceback
from typing import Optional, Tuple, List, Dict, Any, Union
import numpy as np
from osgeo import gdal, ogr, osr
#use rtree for spatial indexing
from rtree.index import Index as RTreeindex
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename
gdal.UseExceptions()
from uraster.classes.sraster import sraster
from uraster.utility import get_polygon_list
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
def run_remap(sFilename_target_mesh,
              sFilename_source_mesh,
              sFilename_source_raster,
              sFilename_raster_mesh,
              dArea_min,
              iFlag_save_clipped_raster_in=0,
              sFolder_raster_out_in=None,
              iFlag_discrete_in=False,
              iFlag_verbose=False,
              iFeature_parallel_threshold=5000):
    """
    Perform zonal statistics by clipping raster data to mesh polygons.

    Main processing method that extracts raster values for each mesh cell polygon
    and computes statistics (mean, min, max, std, sum, count).

    Args:

        sFilename_vector_out (str): Output vector file path with computed statistics
        sFilename_source_mesh_in (str, optional): Input mesh polygon file.
            Defaults to configured target mesh.
        sFilename_source_raster_in (list, optional): List of source raster files.
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



    if iFlag_verbose:
        logger.info("run_remap: Starting input file validation...")
    # check input files

    if os.path.exists(sFilename_source_raster):
        pass
    else:
        logger.error('The raster file does not exist!')
        return

    if iFlag_verbose:
        logger.info(
            f"Checking source mesh file: {os.path.basename(sFilename_source_mesh)}")
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
        # remove the file using the vector driver
        pDriver_vector.DeleteDataSource(sFilename_target_mesh)

    sExtension = os.path.splitext(sFilename_source_raster)[1]
    sName = os.path.basename(sFilename_source_raster)
    sRasterName_no_extension = os.path.splitext(sName)[0]

    if iFlag_verbose:
        logger.info(
            "run_remap: Reading raster metadata and determining processing bounds...")

    # use sraster class to read the raster info
    pRaster = sraster(sFilename_in=sFilename_source_raster)
    pRaster.read_metadata()
    if dPixelWidth is None or pRaster.dResolution_x < dPixelWidth:
        dPixelWidth = pRaster.dResolution_x
    if pPixelHeight is None or abs(pRaster.dResolution_y) < abs(pPixelHeight):
        pPixelHeight = pRaster.dResolution_y
    dMissing_value = pRaster.dNoData


    if iFlag_verbose:
        logger.info("run_remap: Opening mesh dataset and analyzing features...")

    pDateset_source_mesh = pDriver_vector.Open(sFilename_source_mesh, ogr.GA_ReadOnly)
    pLayer_source_mesh = pDateset_source_mesh.GetLayer()
    sProjection_source_wkt = pLayer_source_mesh.GetSpatialRef().ExportToWkt
    #build the rtree index for the polygons for the source mesh
    aPolygon, aArea = get_polygon_list(sFilename_source_mesh,
                                     dArea_min=dArea_min,
                                     iFlag_verbose=iFlag_verbose)
    index_base = RTreeindex()
    for idx, poly in enumerate(aPolygon):
        cellid, wkt = poly
        if wkt is None or wkt == '':
            logger.warning(
                f"run_remap: Warning - Empty geometry for feature ID {cellid}, skipping...")
            continue
        envelope = ogr.CreateGeometryFromWkt(wkt).GetEnvelope()
        left, right, bottom, top = envelope

        # Insert bounding box into spatial index
        # rtree use coordinate order
        pBound = (left, bottom, right, top)
        index_base.insert(cellid, pBound) #can use idx or cellid as the id

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


    pLayer_out.CreateField(ogr.FieldDefn('mean', ogr.OFTReal))
    pLayer_out.CreateField(ogr.FieldDefn('min', ogr.OFTReal))
    pLayer_out.CreateField(ogr.FieldDefn('max', ogr.OFTReal))
    pLayer_out.CreateField(ogr.FieldDefn('std', ogr.OFTReal))


    options = ['COMPRESS=DEFLATE', 'PREDICTOR=2']  # reseverd for future use

    # Pre-compute GDAL options to avoid repeated object creation
    sRemap_method = 'near'
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
        tasks.append((idx, cellid, wkt, sFilename_source_raster,
                     gdal_warp_options_serial, dMissing_value, iFlag_verbose))

    # Choose serial or parallel processing based on threshold
    if n_features <= iFeature_parallel_threshold:
        logger.info(
            f"Feature count ({n_features}) <= threshold ({iFeature_parallel_threshold}); using serial processing")
        for task in tasks:
            feature_idx, cellid, success, payload = _process_task(task)

            if not success:
                failed_features.append(
                    {"feature_id": cellid, "error": payload, "envelope": None})
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
    if iFlag_verbose:
        logger.info(f"Processing completed in {total_time:.2f} seconds")
        logger.info(f"Successfully processed: {successful_features} features")
        logger.info(f"Failed features: {len(failed_features)}")

    if failed_features:
        if iFlag_verbose:
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
            if iFlag_verbose:
                logger.info(f"Failure report saved to: {failure_report_file}")
        except Exception as e:
            logger.error(f"Could not save failure report: {e}")

    return