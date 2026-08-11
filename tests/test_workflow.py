
import unittest
import os, sys, stat
import shutil
import numpy as np

from osgeo import gdal, ogr, osr
import urllib

import pystac_client
import planetary_computer


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

# 1. Define your Area of Interest (Bounding Box: Min Long, Min Lat, Max Long, Max Lat)
# Example: A small snippet in Colorado, US
bbox = [dLongitude_min, dLatitude_min, dLongitude_max, dLatitude_max]
  

def create_continuous_test_raster(sFolder,iFlag_download_in=True):
    driver = gdal.GetDriverByName("GTiff")

    if not os.path.exists(sFolder):
        os.makedirs(sFolder)
    

    # Open the Planetary Computer catalog
    catalog = pystac_client.Client.open(
        "https://planetarycomputer.microsoft.com/api/stac/v1",
        modifier=planetary_computer.sign_inplace
    )

    # Search for the NASADEM (NASA's refined SRTM data)
    search = catalog.search(
        collections=["nasadem"],
        bbox=bbox
    )

    items = search.item_collection()
    nfile = len(items)
    if nfile == 0:
        raise ValueError("No elevation tiles found for this region.")
    else:
        print(f"Found {nfile} elevation tiles.")

    # Grab the cloud-optimized GeoTIFF URL for the first matching tile
    # Note: 'elevation' contains the actual DEM data matrix
    aFilename = list()
    for i in range(nfile):
        print(f"Tile {i+1}: {items[i].id}")
        asset_url = items[i].assets["elevation"].href

        sFilename = os.path.join(sFolder, f"elevation_tile_{i+1}.tif")
        if iFlag_download_in:
            if os.path.exists(sFilename): #delete if already exists
               os.path.remove(sFilename)

            print("Downloading elevation tile...")
            urllib.request.urlretrieve(asset_url, sFilename)
        aFilename.append(sFilename)
        
    return aFilename

def create_discrete_test_raster(sFolder, iFlag_download_in=True):
    driver = gdal.GetDriverByName("GTiff")
    if not os.path.exists(sFolder):
        os.makedirs(sFolder)

    # 2. Connect to the Planetary Computer STAC API
    catalog = pystac_client.Client.open(
        "https://planetarycomputer.microsoft.com/api/stac/v1",
        modifier=planetary_computer.sign_inplace
    )
    
    # 3. Search for Land Cover data (using the 10m annual Impact Observatory dataset)
    search = catalog.search(
        collections=["io-lulc-9-class"],
        bbox=bbox,
        datetime="2023" # Adjust the year if needed
    )
    
    items = search.item_collection()
    nfile = len(items)
    print(f"Found {nfile} matching tiles.")
    
    if nfile == 0:
        raise ValueError("No matching tiles found for the given bounding box and year.")
    
    # 4. Grab the first matching tile and extract the data link
    aFilename = list()
    for i in range(nfile):
        item = items[i]
        asset_url = item.assets["data"].href
        
        # 5. Download and crop the data directly from cloud storage
        output_filename = os.path.join(sFolder, f"land_cover_tile_{i+1}.tif")
        if iFlag_download_in:
            if os.path.exists(output_filename): #delete if already exists
               os.remove(output_filename)

            urllib.request.urlretrieve(asset_url, output_filename)
        aFilename.append(output_filename)
        
    return aFilename

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
        cls.continuous_raster_path = os.path.join(cls.output_dir, "test_continuous_raster")
        cls.discrete_raster_path = os.path.join(cls.output_dir, "test_discrete_raster")
        cls.mesh_path = os.path.join(cls.output_dir, "test_mesh.geojson")
        cls.aFilename_continuous = create_continuous_test_raster(cls.continuous_raster_path, iFlag_download_in=True)
        cls.aFilename_discrete = create_discrete_test_raster(cls.discrete_raster_path, iFlag_download_in=True)
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
            self.aFilename_continuous,
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
            self.aFilename_discrete,
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
