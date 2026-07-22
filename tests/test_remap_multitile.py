import numpy as np
import pytest
from uraster.operation import extract
from uraster import utility


def _area_min(path):
    return utility.rebuild_mesh_topology(
        path, iFlag_verbose_in=False, sField_unique_id="cellid"
    )["area_min"]


def test_two_tiles_are_both_read(make_raster, make_grid_mesh, read_vector, tmp_path):
    # left tile all 4, right tile all 8, equal pixel counts -> mean 6 proves both contribute
    tile_a = make_raster(
        tmp_path / "a.tif", values=np.full((2, 2), 4.0), bbox=(0, 0, 1, 1)
    )
    tile_b = make_raster(
        tmp_path / "b.tif", values=np.full((2, 2), 8.0), bbox=(1, 0, 2, 1)
    )
    mesh = make_grid_mesh(tmp_path / "m.geojson", bbox=(0, 0, 2, 1), nrows=1, ncols=1)
    out = str(tmp_path / "out.geojson")
    extract.run_remap(
        out, mesh, [tile_a, tile_b], _area_min(mesh), iFlag_remap_method_in=1
    )
    rows = read_vector(out)
    assert rows[1]["mean"] == pytest.approx(6.0, abs=1e-6)
