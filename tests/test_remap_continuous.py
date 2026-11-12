import numpy as np
import pytest

from uraster import utility
from uraster.operation import extract


def _area_min(path):
    return utility.rebuild_mesh_topology(
        path, iFlag_verbose_in=False, sField_unique_id="cellid"
    )["area_min"]


def test_constant_raster_mean(make_raster, make_grid_mesh, read_vector, tmp_path):
    raster = make_raster(
        tmp_path / "c.tif", values=np.full((4, 4), 10.0), bbox=(0, 0, 1, 1)
    )
    # single cell matching the raster extent -> captures all 16 pixels
    # (a cell LARGER than the raster zero-fills the outside area and skews the mean)
    mesh = make_grid_mesh(tmp_path / "m.geojson", bbox=(0, 0, 1, 1), nrows=1, ncols=1)
    out = str(tmp_path / "out.geojson")
    extract.run_remap(out, mesh, [raster], _area_min(mesh), iFlag_remap_method_in=1)
    rows = read_vector(out)
    assert rows[1]["mean"] == pytest.approx(10.0, abs=1e-6)


def test_half_and_half_mean(make_raster, make_grid_mesh, read_vector, tmp_path):
    values = np.array([[0, 0, 20, 20]] * 4, dtype=float)  # mean over all 16 == 10
    raster = make_raster(tmp_path / "hh.tif", values=values, bbox=(0, 0, 1, 1))
    mesh = make_grid_mesh(tmp_path / "m.geojson", bbox=(0, 0, 1, 1), nrows=1, ncols=1)
    out = str(tmp_path / "out.geojson")
    extract.run_remap(out, mesh, [raster], _area_min(mesh), iFlag_remap_method_in=1)
    rows = read_vector(out)
    assert rows[1]["mean"] == pytest.approx(10.0, abs=1e-6)


def test_statistics_match_numpy(make_raster, make_grid_mesh, read_vector, tmp_path):
    values = np.arange(16).reshape(4, 4).astype(float)
    raster = make_raster(tmp_path / "s.tif", values=values, bbox=(0, 0, 1, 1))
    mesh = make_grid_mesh(tmp_path / "m.geojson", bbox=(0, 0, 1, 1), nrows=1, ncols=1)
    out = str(tmp_path / "out.geojson")
    extract.run_remap(
        out,
        mesh,
        [raster],
        _area_min(mesh),
        iFlag_remap_method_in=1,
        iFlag_stat_in=True,
    )
    rows = read_vector(out)
    assert rows[1]["mean"] == pytest.approx(float(values.mean()), abs=1e-6)
    assert rows[1]["min"] == pytest.approx(0.0, abs=1e-6)
    assert rows[1]["max"] == pytest.approx(15.0, abs=1e-6)
    assert rows[1]["std"] == pytest.approx(float(values.std()), rel=1e-3)


def test_coarse_raster_to_fine_mesh_all_cells_valid(
    make_raster, make_grid_mesh, read_vector, tmp_path
):
    """
    Test Issue #79: when raster pixels (0.5 deg) are larger than mesh cells (0.1 deg),
    every mesh cell must receive valid statistics (no NaNs) via auto-resampling.
    """
    # 2x2 raster over bbox [0, 1] x [0, 1] -> pixel size = 0.5 deg
    raster = make_raster(
        tmp_path / "coarse.tif", values=np.full((2, 2), 42.0), bbox=(0, 0, 1, 1)
    )
    # 10x10 mesh over bbox [0, 1] x [0, 1] -> cell size = 0.1 deg (smaller than 0.5 deg)
    mesh = make_grid_mesh(
        tmp_path / "fine_mesh.geojson", bbox=(0, 0, 1, 1), nrows=10, ncols=10
    )
    out = str(tmp_path / "out_fine.geojson")
    extract.run_remap(out, mesh, [raster], _area_min(mesh), iFlag_remap_method_in=1)
    rows = read_vector(out)

    assert len(rows) == 100
    for cellid, row in rows.items():
        assert not np.isnan(row["mean"]), f"Cell {cellid} has NaN mean!"
        assert row["mean"] == pytest.approx(42.0, abs=1e-6)


def test_empty_cell_writes_null_not_nan(
    make_raster, make_grid_mesh, read_vector, tmp_path
):
    """Verify that cells with no valid data write explicit NULL (None) rather than NaN."""
    # Raster with all nodata
    raster = make_raster(
        tmp_path / "all_nodata.tif",
        values=np.full((4, 4), -9999.0),
        bbox=(0, 0, 1, 1),
        nodata=-9999.0,
    )
    mesh = make_grid_mesh(tmp_path / "m.geojson", bbox=(0, 0, 1, 1), nrows=1, ncols=1)
    out = str(tmp_path / "out_null.geojson")
    extract.run_remap(out, mesh, [raster], _area_min(mesh), iFlag_remap_method_in=1)
    rows = read_vector(out)
    assert rows[1]["mean"] is None
