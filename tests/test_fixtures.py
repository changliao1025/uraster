import numpy as np
from osgeo import gdal, ogr


def test_make_raster_roundtrip(make_raster, tmp_path):
    path = make_raster(
        tmp_path / "r.tif",
        values=np.arange(16).reshape(4, 4),
        bbox=(0.0, 0.0, 1.0, 1.0),
        epsg=4326,
    )
    ds = gdal.Open(path)
    assert ds.RasterXSize == 4
    assert ds.RasterYSize == 4
    arr = ds.GetRasterBand(1).ReadAsArray()
    # row 0 is the northern edge and holds the first row of values
    assert arr[0, 0] == 0
    assert arr[3, 3] == 15
    ds = None


def test_make_grid_mesh_roundtrip(make_grid_mesh, tmp_path):
    path = make_grid_mesh(tmp_path / "m.geojson", bbox=(0, 0, 3, 3), nrows=3, ncols=3)
    ds = ogr.Open(path)
    layer = ds.GetLayer(0)
    assert layer.GetFeatureCount() == 9
    ids = sorted(f.GetField("cellid") for f in layer)
    assert ids == list(range(1, 10))
    ds = None
