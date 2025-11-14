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
    aPolygon, aArea = get_polygon_list(sFilename_raster_mesh,
                                     iFlag_verbose_in=iFlag_verbose)
    index_mesh = RTreeindex()
    for idx, poly in enumerate(aPolygon):
        cellid, wkt = poly
        if wkt is None or wkt == '':
            logger.warning(
                f"run_remap: Warning - Empty geometry for feature ID {cellid}, skipping...")
            continue
        envelope = ogr.CreateGeometryFromWkt(wkt).GetEnvelope()
        left, right, bottom, top = envelope
        # Insert bounding box into spatial index
        pBound = (left, bottom, right, top)
        index_mesh.insert(idx, pBound) #can use idx or cellid as the id

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



    options = ['COMPRESS=DEFLATE', 'PREDICTOR=2']  # reseverd for future use

    # Pre-compute GDAL options to avoid repeated object creation


    logger.info("run_remap: Starting main feature processing loop...")



    n_features = len(aPolygon)
    max_workers = min(cpu_count(), max(1, n_features))
    logger.info(
        f"Preparing to process {n_features} features (parallel threshold={iFeature_parallel_threshold})")

    start_time = time.time()

    #now we need to find the intersecting polygons between the raster mesh and the source mesh
    for pFeature in pLayer_source_mesh:
        cellid = pFeature.GetFieldAsInteger('cellid')
        pTarget_geometry = pFeature.GetGeometryRef()
        if pTarget_geometry is None:
            logger.warning(
                f"run_remap: Warning - Empty geometry for feature ID {cellid}, skipping...")
            continue
        envelope = pTarget_geometry.GetEnvelope()
        left, right, bottom, top = envelope
        # Query spatial index for candidate intersecting polygons
        candidate_idxs = list(index_mesh.intersection((left, bottom, right, top)))
        # Further process candidates to find actual intersections
        for idx in candidate_idxs:
            raster_cellid, raster_wkt = aPolygon[idx]
            raster_geometry = ogr.CreateGeometryFromWkt(raster_wkt)
            #first check whether the mesh is inside the target polygon
            if pTarget_geometry.Contains(raster_geometry):
                # keep the raster geometry for further processing
                pass
            else:
                if pTarget_geometry.Intersects(raster_geometry): #both intersect and touching
                    # get the intersected geometry
                    pIntersected_geometry = pTarget_geometry.Intersection(
                        raster_geometry)
                    #should be a polygon geometry?
                    sGeometryName = pIntersected_geometry.GetGeometryName()
                    if sGeometryName == 'POLYGON':
                        pass
                else:
                    continue  # no intersection, skip





    # flush and close output
    pDataset_out.FlushCache()
    pDataset_out = None

    # Clean up spatial reference objects to prevent memory leaks
    pSpatialRef_target = None

    # Report processing summary
    total_time = time.time() - start_time
    if iFlag_verbose:
        logger.info(f"Processing completed in {total_time:.2f} seconds")




    return