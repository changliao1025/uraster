from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import Pool, cpu_count
import os
import time
from typing import Optional, Tuple, List, Dict, Any, Union
import numpy as np
from osgeo import gdal, ogr, osr

# use rtree for spatial indexing
from rtree.index import Index as RTreeindex
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area

gdal.UseExceptions()
from uraster.classes.sraster import sraster
from uraster.utility import get_polygon_list

# Try to import psutil for memory monitoring (optional)
try:
    import psutil

    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False

from uraster.utility import setup_logger

logger = setup_logger(__name__.split(".")[-1])
crs = "EPSG:4326"

# Initialize GDAL drivers with error handling
try:
    pDriver_geojson = ogr.GetDriverByName("GeoJSON")
    pDriver_shp = ogr.GetDriverByName("ESRI Shapefile")
    if pDriver_geojson is None or pDriver_shp is None:
        raise RuntimeError("Failed to initialize required GDAL drivers")
except Exception as e:
    logger.error(f"Error initializing GDAL drivers: {e}")
    raise

# Constants for processing thresholds
IDL_LONGITUDE_THRESHOLD = 100  # Degrees - threshold for detecting IDL crossing
WARP_TIMEOUT_SECONDS = 30  # Seconds - timeout for GDAL Warp operations
PROGRESS_REPORT_INTERVAL = 5  # Report progress every N features
MAX_CONSECUTIVE_FAILURES = 10  # Maximum consecutive failures before stopping
HEARTBEAT_INTERVAL = 5  # Seconds between heartbeat logs during long operations


def run_remap(
    sFilename_target_mesh,
    sFilename_source_mesh,
    sFilename_source_raster,
    sFilename_raster_mesh,
    iFlag_save_clipped_raster_in=0,
    sFolder_raster_out_in=None,
    iFlag_discrete_in=False,
    iFlag_verbose_in=False,
    iFeature_parallel_threshold=5000,
    sField_unique_id="cellid",
):
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
        iFlag_verbose_in (bool, optional): If True, print detailed progress messages.
            If False, only print error messages. Default is False.

    Returns:
        None

    Note:
        - Handles IDL-crossing polygons automatically
        - Generates failure report for problematic features
        - Supports multiple input rasters (uses first in list)
    """

    if iFlag_verbose_in:
        logger.info("run_remap: Starting input file validation...")
    # check input files

    if os.path.exists(sFilename_source_raster):
        pass
    else:
        logger.error("The raster file does not exist!")
        return

    if iFlag_verbose_in:
        logger.info(
            f"Checking source mesh file: {os.path.basename(sFilename_source_mesh)}"
        )
    if os.path.exists(sFilename_source_mesh):
        pass
    else:
        logger.error("The vector mesh file does not exist!")
        return
    if iFlag_verbose_in:
        logger.info("Input file validation completed successfully")

    # Determine output vector format from filename extension
    pDriver_vector = get_vector_driver_from_filename(sFilename_target_mesh)

    if os.path.exists(sFilename_target_mesh):
        # remove the file using the vector driver
        pDriver_vector.DeleteDataSource(sFilename_target_mesh)

    sExtension = os.path.splitext(sFilename_source_raster)[1].lstrip(".")
    sName = os.path.basename(sFilename_source_raster)
    sRasterName_no_extension = os.path.splitext(sName)[0]

    if iFlag_verbose_in:
        logger.info(
            "run_remap: Reading raster metadata and determining processing bounds..."
        )

    # use sraster class to read the raster info
    pRaster = sraster(sFilename_in=sFilename_source_raster)
    pRaster.read_metadata()
    pRaster_data = pRaster.read_data(iBand=1)
    # Initialize pixel resolution variables
    dPixelWidth = pRaster.dResolution_x
    pPixelHeight = pRaster.dResolution_y
    dMissing_value = pRaster.dNoData

    if iFlag_verbose_in:
        logger.info("run_remap: Opening mesh dataset and analyzing features...")

    pDateset_source_mesh = pDriver_vector.Open(sFilename_source_mesh, 0)
    pLayer_source_mesh = pDateset_source_mesh.GetLayer()
    sProjection_source_wkt = pLayer_source_mesh.GetSpatialRef().ExportToWkt
    # build the rtree index for the polygons for the source mesh
    aPolygon, aArea, sProjection_source_wkt = get_polygon_list(
        sFilename_raster_mesh,
        iFlag_verbose_in=iFlag_verbose_in,
        sField_unique_id=sField_unique_id,
    )
    index_raster_mesh = RTreeindex()  # build the spatial index for the raster mesh cell
    for idx, poly in enumerate(aPolygon):
        cellid, wkt = poly
        if wkt is None or wkt == "":
            logger.warning(
                f"run_remap: Warning - Empty geometry for feature ID {cellid}, skipping..."
            )
            continue
        envelope = ogr.CreateGeometryFromWkt(wkt).GetEnvelope()
        left, right, bottom, top = envelope
        # Insert bounding box into spatial index
        pBound = (left, bottom, right, top)
        index_raster_mesh.insert(idx, pBound)  # can use idx or cellid as the id

    pSpatialRef_target = osr.SpatialReference()
    pSpatialRef_target.ImportFromWkt(sProjection_source_wkt)

    # create a polygon feature to save the output
    pDataset_out = pDriver_vector.CreateDataSource(sFilename_target_mesh)
    pLayer_out = pDataset_out.CreateLayer("uraster", pSpatialRef_target, ogr.wkbPolygon)
    pLayer_defn_out = pLayer_out.GetLayerDefn()
    pFeature_out = ogr.Feature(pLayer_defn_out)
    pLayer_out.CreateField(ogr.FieldDefn(sField_unique_id, ogr.OFTInteger))
    # define a field
    pField = ogr.FieldDefn("area", ogr.OFTReal)
    pField.SetWidth(32)
    pField.SetPrecision(2)
    pLayer_out.CreateField(pField)
    # in the future, we will also copy other attributes from the input geojson file
    pLayer_out.CreateField(ogr.FieldDefn("mean", ogr.OFTReal))
    options = ["COMPRESS=DEFLATE", "PREDICTOR=2"]  # reseverd for future use
    # Pre-compute GDAL options to avoid repeated object creation
    logger.info("run_remap: Starting main feature processing loop...")
    n_raster_features = len(aPolygon)
    max_workers = min(cpu_count(), max(1, n_raster_features))
    n_source_features = pLayer_source_mesh.GetFeatureCount()
    logger.info(f"Total number of source mesh features to process: {n_source_features}")
    start_time = time.time()
    # now we need to find the intersecting polygons between the raster mesh and the source mesh
    for pFeature in pLayer_source_mesh:
        cellid = pFeature.GetFieldAsInteger(sField_unique_id)
        pTarget_geometry = pFeature.GetGeometryRef()
        # get name
        if pTarget_geometry is None:
            logger.warning(
                f"run_remap: Warning - Empty geometry for feature ID {cellid}, skipping..."
            )
            continue
        aCoords = get_geometry_coordinates(pTarget_geometry)
        sGeometryName = pTarget_geometry.GetGeometryName()
        aMesh_cell_within = list()
        aArea_ratio = list()
        aMesh_cell_intersect = list()
        if sGeometryName == "POLYGON":
            dArea_total_source = calculate_polygon_area(aCoords[:, 0], aCoords[:, 1])
            envelope = pTarget_geometry.GetEnvelope()
            left, right, bottom, top = envelope
            # Query spatial index for candidate intersecting polygons
            candidate_idxs = list(
                index_raster_mesh.intersection((left, bottom, right, top))
            )
            # Further process candidates to find actual intersections

            for idx in candidate_idxs:
                id, raster_wkt = aPolygon[idx]
                raster_geometry = ogr.CreateGeometryFromWkt(raster_wkt)
                # first check whether the mesh is inside the target polygon
                if pTarget_geometry.Contains(raster_geometry):
                    # keep the raster geometry for further processing
                    aMesh_cell_within.append(idx)
                    pass
                else:
                    if pTarget_geometry.Intersects(
                        raster_geometry
                    ):  # both intersect and touching
                        # get the intersected geometry
                        pIntersected_geometry = pTarget_geometry.Intersection(
                            raster_geometry
                        )
                        # should be a polygon geometry?
                        sGeometryName = pIntersected_geometry.GetGeometryName()
                        if sGeometryName == "POLYGON":
                            # Get the area of the intersected polygon
                            dArea_raster = aArea[idx]
                            aCoords_intersect = get_geometry_coordinates(
                                pIntersected_geometry
                            )
                            dArea_intersect = calculate_polygon_area(
                                aCoords_intersect[:, 0], aCoords_intersect[:, 1]
                            )
                            # check the area ratio
                            aMesh_cell_intersect.append(idx)
                            aArea_ratio.append(dArea_intersect / dArea_raster)
                            pass
                    else:
                        continue  # no intersection, skip

        elif sGeometryName == "MULTIPOLYGON":
            dArea_total_source = 0.0
            for i in range(pTarget_geometry.GetGeometryCount()):
                pSub_geometry = pTarget_geometry.GetGeometryRef(i)
                # check validity
                if not pSub_geometry.IsValid():
                    # flatten to 2d
                    pSub_geometry.FlattenTo2D()
                    logger.warning(
                        f"Problematic geometry at coordinates: {pSub_geometry.ExportToWkt()}"
                    )
                    logger.warning("Invalid geometry detected, attempting to fix...")
                    continue

                aCoords_sub = get_geometry_coordinates(pSub_geometry)
                dArea_sub = calculate_polygon_area(aCoords_sub[:, 0], aCoords_sub[:, 1])
                dArea_total_source += dArea_sub

                # now we need to find the intersecting raster mesh cells for each sub-geometry
                envelope = pSub_geometry.GetEnvelope()
                left, right, bottom, top = envelope
                # Query spatial index for candidate intersecting polygons
                candidate_idxs = list(
                    index_raster_mesh.intersection((left, bottom, right, top))
                )
                # Further process candidates to find actual intersections
                for idx in candidate_idxs:
                    id, raster_wkt = aPolygon[idx]
                    raster_geometry = ogr.CreateGeometryFromWkt(raster_wkt)
                    # first check whether the mesh is inside the target polygon
                    if pSub_geometry.Contains(raster_geometry):
                        # keep the raster geometry for further processing
                        if idx not in aMesh_cell_within:
                            aMesh_cell_within.append(idx)
                        pass
                    else:
                        if pSub_geometry.Intersects(
                            raster_geometry
                        ):  # both intersect and touching
                            # get the intersected geometry
                            pIntersected_geometry = pSub_geometry.Intersection(
                                raster_geometry
                            )
                            # should be a polygon geometry?
                            sGeometryName = pIntersected_geometry.GetGeometryName()
                            if sGeometryName == "POLYGON":
                                # Get the area of the intersected polygon
                                dArea_raster = aArea[idx]
                                aCoords_intersect = get_geometry_coordinates(
                                    pIntersected_geometry
                                )
                                dArea_intersect = calculate_polygon_area(
                                    aCoords_intersect[:, 0], aCoords_intersect[:, 1]
                                )
                                # check the area ratio
                                aMesh_cell_intersect.append(idx)
                                aArea_ratio.append(dArea_intersect / dArea_raster)
                                pass
                        else:
                            continue  # no intersection, skip

        # create the output feature
        pFeature_out = ogr.Feature(pLayer_defn_out)
        # set the id field
        pFeature_out.SetField(sField_unique_id, cellid)
        # set the geometry as the target geometry
        pFeature_out.SetGeometry(pTarget_geometry.Clone())

        # now we can calculate the weighted area for the target

        dArea_total = dArea_total_source
        dArea_check = 0.0
        # the weighted mean equation is:
        # $$
        # \text{Weighted Mean} = \frac{\sum_{i=1}^{N} (\text{Value}_i \times \text{Area}_i)}{\text{Total Area}}
        # $$
        dWeighted_sum = 0.0
        # first process the within cells
        for idx in aMesh_cell_within:
            # convert the ids to row and column indices
            nCol = idx % pRaster.ncolumn
            nRow = idx // pRaster.ncolumn
            # the idx start from the lower left corner, so we need to convert it to the upper left corner
            nRow_converted = pRaster.nrow - 1 - nRow
            # get the raster value
            dRaster_value = pRaster_data[nRow_converted, nCol]
            dArea_mesh_cell = aArea[idx]
            if dRaster_value == dMissing_value:
                continue
            else:
                dWeighted_sum += dRaster_value * dArea_mesh_cell
                dArea_check += dArea_mesh_cell

        # then process the intersected cells
        for k, idx in enumerate(aMesh_cell_intersect):
            # convert the ids to row and column indices
            nCol = idx % pRaster.ncolumn
            nRow = idx // pRaster.ncolumn
            # the idx start from the lower left corner, so we need to convert it to the upper left corner
            nRow_converted = pRaster.nrow - 1 - nRow
            # get the raster value
            dRaster_value = pRaster_data[nRow_converted, nCol]
            dArea_mesh_cell = aArea[idx]
            dArea_ratio = aArea_ratio[k]
            dArea_intersected = dArea_mesh_cell * dArea_ratio
            if dRaster_value == dMissing_value:
                continue
            else:
                dWeighted_sum += dRaster_value * dArea_intersected
                dArea_check += dArea_intersected

        if dArea_check == 0:
            dWeighted_mean = None
        else:
            # the difference between dArea_check and dArea_total should be small
            if abs(dArea_check - dArea_total) / dArea_total > 0.01:
                logger.warning(
                    f"run_remap: Warning - Area check failed for feature ID {cellid}..."
                )
                if dArea_check == 0:
                    dWeighted_mean = None
                else:
                    dWeighted_mean = dWeighted_sum / dArea_check
            else:
                dWeighted_mean = dWeighted_sum / dArea_total

        # set the area field
        pFeature_out.SetField("area", dArea_total_source)
        # set the mean field
        if dWeighted_mean is not None:
            pFeature_out.SetField("mean", dWeighted_mean)
        else:
            pFeature_out.SetField("mean", None)
        # create the feature in the output layer
        pLayer_out.CreateFeature(pFeature_out)
        pFeature_out = None  # destroy the feature to free memory

    # flush and close output
    pDataset_out.FlushCache()
    pDataset_out = None

    # Clean up spatial reference objects to prevent memory leaks
    pSpatialRef_target = None

    # Report processing summary
    total_time = time.time() - start_time
    if iFlag_verbose_in:
        logger.info(f"Processing completed in {total_time:.2f} seconds")

    return
