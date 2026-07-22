import os
import numpy as np
import pytest
from osgeo import gdal, ogr, osr

gdal.UseExceptions()


def _make_raster(path, values, bbox, epsg=4326, nodata=None, dtype=gdal.GDT_Float32):
    values = np.asarray(values, dtype=float)
    nrow, ncol = values.shape
    min_x, min_y, max_x, max_y = bbox
    res_x = (max_x - min_x) / ncol
    res_y = (max_y - min_y) / nrow
    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create(str(path), ncol, nrow, 1, dtype)
    # top-left x, w-e res, 0, top-left y, 0, n-s res (negative)
    ds.SetGeoTransform((min_x, res_x, 0.0, max_y, 0.0, -res_y))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    ds.SetProjection(srs.ExportToWkt())
    band = ds.GetRasterBand(1)
    if nodata is not None:
        band.SetNoDataValue(float(nodata))
    band.WriteArray(values)
    band.FlushCache()
    ds = None
    return str(path)


def _make_grid_mesh(path, bbox, nrows, ncols, epsg=4326, id_field="cellid"):
    min_x, min_y, max_x, max_y = bbox
    dx = (max_x - min_x) / ncols
    dy = (max_y - min_y) / nrows
    driver = ogr.GetDriverByName("GeoJSON")
    if os.path.exists(path):
        driver.DeleteDataSource(str(path))
    ds = driver.CreateDataSource(str(path))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    layer = ds.CreateLayer("mesh", srs, ogr.wkbPolygon)
    layer.CreateField(ogr.FieldDefn(id_field, ogr.OFTInteger))
    cid = 1
    for i in range(nrows):
        for j in range(ncols):
            x0 = min_x + j * dx
            x1 = x0 + dx
            y0 = min_y + i * dy
            y1 = y0 + dy
            ring = ogr.Geometry(ogr.wkbLinearRing)
            for (x, y) in [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]:
                ring.AddPoint(x, y)
            ring.CloseRings()
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)
            feat = ogr.Feature(layer.GetLayerDefn())
            feat.SetField(id_field, cid)
            feat.SetGeometry(poly)
            layer.CreateFeature(feat)
            feat = None
            cid += 1
    ds = None
    return str(path)


def _make_polygon_mesh(path, polygons, epsg=4326, id_field="cellid"):
    driver = ogr.GetDriverByName("GeoJSON")
    if os.path.exists(path):
        driver.DeleteDataSource(str(path))
    ds = driver.CreateDataSource(str(path))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    layer = ds.CreateLayer("mesh", srs, ogr.wkbPolygon)
    layer.CreateField(ogr.FieldDefn(id_field, ogr.OFTInteger))
    for idx, coords in enumerate(polygons, start=1):
        ring = ogr.Geometry(ogr.wkbLinearRing)
        for (x, y) in coords:
            ring.AddPoint(float(x), float(y))
        ring.CloseRings()
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)
        feat = ogr.Feature(layer.GetLayerDefn())
        feat.SetField(id_field, idx)
        feat.SetGeometry(poly)
        layer.CreateFeature(feat)
        feat = None
    ds = None
    return str(path)


def _read_vector(path):
    ds = ogr.Open(str(path))
    assert ds is not None, f"could not open {path}"
    layer = ds.GetLayer(0)
    defn = layer.GetLayerDefn()
    names = [defn.GetFieldDefn(i).GetName() for i in range(defn.GetFieldCount())]
    rows = {}
    for feat in layer:
        rec = {n: feat.GetField(n) for n in names}
        geom = feat.GetGeometryRef()
        rec["__geom__"] = geom.Clone() if geom is not None else None
        rows[rec.get("cellid")] = rec
    ds = None
    return rows


@pytest.fixture
def make_raster():
    return _make_raster


@pytest.fixture
def make_grid_mesh():
    return _make_grid_mesh


@pytest.fixture
def make_polygon_mesh():
    return _make_polygon_mesh


@pytest.fixture
def read_vector():
    return _read_vector
