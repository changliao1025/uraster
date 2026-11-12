import numpy as np
import pytest

from uraster.classes.sraster import SRaster
from uraster.operation import intersect


def _raster_mesh(raster_path):
    r = SRaster(sFilename_in=raster_path)
    r.read_metadata()
    return r.create_raster_mesh()


def test_weighted_constant(make_raster, make_grid_mesh, read_vector, tmp_path):
    raster = make_raster(
        tmp_path / "wc.tif", values=np.full((4, 4), 10.0), bbox=(0, 0, 1, 1)
    )
    mesh = make_grid_mesh(
        tmp_path / "m.geojson", bbox=(0.1, 0.1, 0.9, 0.9), nrows=1, ncols=1
    )
    out = str(tmp_path / "out.geojson")
    intersect.run_remap(out, mesh, raster, _raster_mesh(raster), iFlag_verbose_in=False)
    rows = read_vector(out)
    # rel tol (not 1e-6): float32 area-weighting over large area values rounds ~1e-4
    assert rows[1]["mean"] == pytest.approx(10.0, rel=1e-3)


def test_weighted_half_and_half_near_equator(
    make_raster, make_grid_mesh, read_vector, tmp_path
):
    # near-equator, tiny extent -> pixels are ~equal area, so area-weighted mean ~ simple mean
    values = np.array([[0, 0, 20, 20]] * 4, dtype=float)
    raster = make_raster(tmp_path / "wh.tif", values=values, bbox=(0, 0, 0.01, 0.01))
    mesh = make_grid_mesh(
        tmp_path / "m.geojson", bbox=(0.001, 0.001, 0.009, 0.009), nrows=1, ncols=1
    )
    out = str(tmp_path / "out.geojson")
    intersect.run_remap(out, mesh, raster, _raster_mesh(raster), iFlag_verbose_in=False)
    rows = read_vector(out)
    assert rows[1]["mean"] == pytest.approx(10.0, rel=1e-3)


def test_weighted_nan_nodata_not_contaminating(
    make_raster, make_grid_mesh, read_vector, tmp_path
):
    """Test Issue #80 Point 4: NaN nodata pixels should be skipped and not produce NaN cells."""
    values = np.array([[10.0, 10.0], [np.nan, 10.0]], dtype=np.float32)  # one NaN pixel
    raster = make_raster(
        tmp_path / "nan_nodata.tif", values=values, bbox=(0, 0, 1, 1), nodata=np.nan
    )
    mesh = make_grid_mesh(
        tmp_path / "m_nan.geojson", bbox=(0.1, 0.1, 0.9, 0.9), nrows=1, ncols=1
    )
    out = str(tmp_path / "out_nan.geojson")
    intersect.run_remap(out, mesh, raster, _raster_mesh(raster), iFlag_verbose_in=False)
    rows = read_vector(out)
    assert rows[1]["mean"] is not None
    assert not np.isnan(rows[1]["mean"])
    # All valid pixels are 10.0, so valid-area mean must be 10.0
    assert rows[1]["mean"] == pytest.approx(10.0, rel=1e-3)


def test_weighted_valid_area_denominator_no_bias(
    make_raster, make_grid_mesh, read_vector, tmp_path
):
    """Test Issue #80 Point 5: mean divides by valid area (dArea_check), so partial missing is not biased low."""
    # 10x10 raster, all 20.0 except one pixel missing
    values = np.full((10, 10), 20.0, dtype=np.float32)
    values[0, 0] = -9999.0
    raster = make_raster(
        tmp_path / "partial_missing.tif",
        values=values,
        bbox=(0, 0, 1, 1),
        nodata=-9999.0,
    )
    mesh = make_grid_mesh(
        tmp_path / "m_partial.geojson", bbox=(0.0, 0.0, 1.0, 1.0), nrows=1, ncols=1
    )
    out = str(tmp_path / "out_partial.geojson")
    intersect.run_remap(out, mesh, raster, _raster_mesh(raster), iFlag_verbose_in=False)
    rows = read_vector(out)
    # The valid pixels are all 20.0, so the mean over valid area must be exactly 20.0 (not biased low)
    assert rows[1]["mean"] == pytest.approx(20.0, rel=1e-3)


def test_weighted_different_source_and_target_formats(
    make_raster, make_grid_mesh, read_vector, tmp_path
):
    """Test Issue #80 Point 3: GeoJSON source mesh with GPKG target mesh should not fail driver opening."""
    from osgeo import ogr

    target_ext = "gpkg" if ogr.GetDriverByName("GPKG") is not None else "geojson"
    raster = make_raster(
        tmp_path / "diff_fmt.tif", values=np.full((4, 4), 10.0), bbox=(0, 0, 1, 1)
    )
    mesh = make_grid_mesh(
        tmp_path / "m_src.geojson", bbox=(0.1, 0.1, 0.9, 0.9), nrows=1, ncols=1
    )
    out = str(tmp_path / f"out_tgt.{target_ext}")
    intersect.run_remap(out, mesh, raster, _raster_mesh(raster), iFlag_verbose_in=False)
    rows = read_vector(out)
    assert rows[1]["mean"] == pytest.approx(10.0, rel=1e-3)


def test_sraster_mesh_driver_availability(make_raster, tmp_path):
    """Test Issue #80 Point 1: sraster.sFilename_mesh selects a driver available in current GDAL."""
    from osgeo import ogr

    raster = make_raster(
        tmp_path / "drv_test.tif", values=np.full((2, 2), 1.0), bbox=(0, 0, 1, 1)
    )
    r = SRaster(sFilename_in=raster)
    mesh_path = r.sFilename_mesh
    ext = mesh_path.split(".")[-1]
    driver_map = {
        "geoparquet": "Parquet",
        "fgb": "FlatGeobuf",
        "gpkg": "GPKG",
        "geojson": "GeoJSON",
    }
    expected_driver_name = driver_map[ext]
    assert ogr.GetDriverByName(expected_driver_name) is not None


def test_sraster_mesh_driver_fallback(monkeypatch, make_raster, tmp_path):
    """Test dynamic fallback of mesh vector drivers based on available OGR drivers."""
    raster = make_raster(
        tmp_path / "fallback_test.tif", values=np.full((2, 2), 1.0), bbox=(0, 0, 1, 1)
    )

    # 1. Parquet available -> .geoparquet
    monkeypatch.setattr("uraster.classes.sraster.has_parquet_support", lambda: True)
    r_pq = SRaster(sFilename_in=raster)
    assert r_pq.sFilename_mesh.endswith("_mesh.geoparquet")

    # 2. No Parquet, but FlatGeobuf available -> .fgb
    monkeypatch.setattr("uraster.classes.sraster.has_parquet_support", lambda: False)
    monkeypatch.setattr(
        "uraster.classes.sraster.get_available_vector_formats",
        lambda: {".fgb": "FlatGeobuf", ".gpkg": "GPKG", ".geojson": "GeoJSON"},
    )
    r_fgb = SRaster(sFilename_in=raster)
    assert r_fgb.sFilename_mesh.endswith("_mesh.fgb")

    # 3. Only GPKG and GeoJSON available -> .gpkg
    monkeypatch.setattr(
        "uraster.classes.sraster.get_available_vector_formats",
        lambda: {".gpkg": "GPKG", ".geojson": "GeoJSON"},
    )
    r_gpkg = SRaster(sFilename_in=raster)
    assert r_gpkg.sFilename_mesh.endswith("_mesh.gpkg")

    # 4. Only GeoJSON available -> .geojson
    monkeypatch.setattr(
        "uraster.classes.sraster.get_available_vector_formats",
        lambda: {".geojson": "GeoJSON"},
    )
    r_geojson = SRaster(sFilename_in=raster)
    assert r_geojson.sFilename_mesh.endswith("_mesh.geojson")


def test_sraster_raster_driver_metadata(make_raster, tmp_path):
    """Test that sraster correctly detects raster driver and format."""
    raster = make_raster(
        tmp_path / "meta_test.tif", values=np.full((2, 2), 1.0), bbox=(0, 0, 1, 1)
    )
    r = SRaster(sFilename_in=raster)
    r.read_metadata()
    assert r.sDriver_raster == "GTiff"
    assert r.sDriver == "GTiff"
