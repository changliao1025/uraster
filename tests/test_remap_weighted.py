import numpy as np
import pytest
from uraster.operation import intersect
from uraster.classes.sraster import sraster


def _raster_mesh(raster_path):
    r = sraster(sFilename_in=raster_path)
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
