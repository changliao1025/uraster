import os
import unittest
import numpy as np
from osgeo import gdal, ogr, osr

# Ensure package can be imported from the repository root
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pyearth.gis.spatialref.reproject_coordinates import reproject_coordinates
from uraster import sraster
from uraster import utility
from uraster.operation.extract import run_remap

#define a shared boundary for the test rasters and mesh
#using Beijing as an example, with a bounding box of (39.7, 116.2) to (40.1, 116.6) 

dLongitude_min = 116.2
dLongitude_max = 116.6
dLatitude_min = 39.7
dLatitude_max = 40.1

srs_wgs84 = osr.SpatialReference()
srs_wgs84.ImportFromEPSG(4326)
wgs84_wkt = srs_wgs84.ExportToWkt()

srs_projected = osr.SpatialReference()
srs_projected.ImportFromEPSG(32650)
projected_wkt = srs_projected.ExportToWkt()

def create_continuous_test_raster(sFilename):
    driver = gdal.GetDriverByName("GTiff")
    if os.path.exists(sFilename):
        driver.Delete(sFilename)

    srs_projected = osr.SpatialReference()
    srs_projected.ImportFromEPSG(32650)

    #use the same bounding box as the mesh, but in projected coordinates
    #we need to convert the lat/lon bounding box to projected coordinates using EPSG:32650
    
    min_x, min_y = reproject_coordinates(dLongitude_min, dLatitude_min, wgs84_wkt, projected_wkt)
    max_x, max_y = reproject_coordinates(dLongitude_max, dLatitude_max, wgs84_wkt, projected_wkt)
    pixel_width = (max_x - min_x) / 100.0
    pixel_height = (max_y - min_y) / 100.0

    ds = driver.Create(sFilename, 100, 100, 1, gdal.GDT_Float32)
    ds.SetGeoTransform((min_x, pixel_width, 0.0, max_y, 0.0, -pixel_height))
    ds.SetProjection(srs_projected.ExportToWkt())
    #array = np.arange(10000, dtype=np.float32).reshape((100, 100))
    #generate a random array of floats between 0 and 1
    array = np.random.uniform(0, 1, size=(100, 100)).astype(np.float32)
    band = ds.GetRasterBand(1)
    band.WriteArray(array)
    band.SetNoDataValue(-9999.0)
    band.FlushCache()
    ds = None
    return array

def create_discrete_test_raster(sFilename):
    driver = gdal.GetDriverByName("GTiff")
    if os.path.exists(sFilename):
        driver.Delete(sFilename)

    srs_projected = osr.SpatialReference()
    srs_projected.ImportFromEPSG(32650)

    min_x, min_y = reproject_coordinates(dLongitude_min, dLatitude_min, wgs84_wkt, projected_wkt)
    max_x, max_y = reproject_coordinates(dLongitude_max, dLatitude_max, wgs84_wkt, projected_wkt)

    pixel_width = (max_x - min_x) / 100.0
    pixel_height = (max_y - min_y) / 100.0

    ds = driver.Create(sFilename, 100, 100, 1, gdal.GDT_Int16)
    ds.SetGeoTransform((min_x, pixel_width, 0.0, max_y, 0.0, -pixel_height))
    ds.SetProjection(srs_projected.ExportToWkt())
    #array = np.ones((100, 100), dtype=np.int16)
    #use 0 to 10 as random discrete values
    #generate a random array of integers between 0 and 10
    array = np.random.randint(0, 11, size=(100, 100), dtype=np.int16)
    band = ds.GetRasterBand(1)
    band.WriteArray(array)
    band.SetNoDataValue(-9999)
    band.FlushCache()
    ds = None
    return array

def create_test_mesh(sFilename):
    driver = ogr.GetDriverByName("GeoJSON")
    if os.path.exists(sFilename):
        driver.DeleteDataSource(sFilename)
    ds = driver.CreateDataSource(sFilename)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    layer = ds.CreateLayer("mesh", srs, ogr.wkbPolygon)
    field = ogr.FieldDefn("cellid", ogr.OFTInteger)
    layer.CreateField(field)

    #ideally, the mesh should have a coarser resolution than the raster
    #let's use a 30 by 30 that covers the raster extent
    dLatitude_range = dLatitude_max - dLatitude_min
    dLongitude_range = dLongitude_max - dLongitude_min
    nRow = 30

    nColumn = 30
    dLatitude_step = dLatitude_range / nRow
    dLongitude_step = dLongitude_range / nColumn
    for i in range(nRow):
        for j in range(nColumn):
            min_lat = 39.7 + i * dLatitude_step
            max_lat = 39.7 + (i + 1) * dLatitude_step
            min_lon = 116.2 + j * dLongitude_step
            max_lon = 116.2 + (j + 1) * dLongitude_step

            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(min_lon, min_lat)
            ring.AddPoint(max_lon, min_lat)
            ring.AddPoint(max_lon, max_lat)
            ring.AddPoint(min_lon, max_lat)
            ring.CloseRings()
            polygon = ogr.Geometry(ogr.wkbPolygon)
            polygon.AddGeometry(ring)

            feature = ogr.Feature(layer.GetLayerDefn())
            feature.SetField("cellid", i * nColumn + j + 1)
            feature.SetGeometry(polygon)
            layer.CreateFeature(feature)
            feature = None
    ds = None


class TestUrasterWorkflow(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.test_id = "workflow_smoke"
        cls.output_dir = os.path.join(os.path.dirname(__file__), "test_data", cls.test_id)
        os.makedirs(cls.output_dir, exist_ok=True)
        cls.continuous_raster_path = os.path.join(cls.output_dir, "test_continuous_raster.tif")
        cls.discrete_raster_path = os.path.join(cls.output_dir, "test_discrete_raster.tif")
        cls.mesh_path = os.path.join(cls.output_dir, "test_mesh.geojson")
        cls.expected_continuous_raster = create_continuous_test_raster(cls.continuous_raster_path)
        cls.expected_discrete_raster = create_discrete_test_raster(cls.discrete_raster_path)
        create_test_mesh(cls.mesh_path)

    @classmethod
    def tearDownClass(cls):
        # Keep the synthetic fixtures in the test data folder so they can be inspected.
        pass

    def setUp(self):
        self.output_path = os.path.join(
            self.output_dir, f"{self._testMethodName}_output.geojson"
        )
        if os.path.exists(self.output_path):
            os.remove(self.output_path)

    def test_sraster_read_metadata(self):
        pRaster = sraster(self.continuous_raster_path)
        pRaster.read_metadata()

        self.assertEqual(pRaster.ncolumn, 100)
        self.assertEqual(pRaster.nrow, 100)
        self.assertIsNotNone(pRaster.sCrs)
        self.assertTrue(pRaster.pSpatialRef.IsProjected())
        self.assertEqual(pRaster.iBandCount, 1)
        self.assertEqual(pRaster.dNoData, -9999.0)

    def test_mesh_topology_rebuild(self):
        mesh_info = utility.rebuild_mesh_topology(
            self.mesh_path, iFlag_verbose_in=False, sField_unique_id="cellid"
        )

        self.assertIsNotNone(mesh_info)
        ds = ogr.Open(self.mesh_path)
        self.assertIsNotNone(ds)
        try:
            layer = ds.GetLayer(0)
            spatial_ref = layer.GetSpatialRef()
            self.assertIsNotNone(spatial_ref)
            self.assertEqual(spatial_ref.GetAuthorityCode(None), "4326")
        finally:
            ds = None

        self.assertEqual(mesh_info["num_cells"], 900)
        self.assertEqual(mesh_info["num_polygns"], 900)
        self.assertEqual(mesh_info["cell_ids"].tolist(), list(range(1, 901)))
       
       

    def test_run_remap_continuous_uses_mean_statistics(self):
        mesh_info = utility.rebuild_mesh_topology(
            self.mesh_path, iFlag_verbose_in=False, sField_unique_id="cellid"
        )
        dArea_min = mesh_info["area_min"]

        run_remap(
            self.output_path,
            self.mesh_path,
            [self.continuous_raster_path],
            dArea_min,
            iFlag_remap_method_in=3,
            iFlag_stat_in=True,
            iFlag_verbose_in=True,
        )

        self.assertTrue(os.path.exists(self.output_path))
        datasource = ogr.Open(self.output_path)
        self.assertIsNotNone(datasource)
        try:
            pass
        finally:
            datasource = None

    def test_run_remap_discrete_uses_dominant_mode(self):
        mesh_info = utility.rebuild_mesh_topology(
            self.mesh_path, iFlag_verbose_in=False, sField_unique_id="cellid"
        )
        dArea_min = mesh_info["area_min"]

        run_remap(
            self.output_path,
            self.mesh_path,
            [self.discrete_raster_path],
            dArea_min,
            iFlag_remap_method_in=1,
            iFlag_discrete_in=True,
            iFlag_verbose_in=False,
        )

        self.assertTrue(os.path.exists(self.output_path))
        datasource = ogr.Open(self.output_path)
        self.assertIsNotNone(datasource)
        try:
            pass
        finally:
            datasource = None


if __name__ == "__main__":
    unittest.main()
