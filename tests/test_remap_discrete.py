import numpy as np
import pytest
from uraster.operation import extract
from uraster import utility


def _area_min(path):
    return utility.rebuild_mesh_topology(
        path, iFlag_verbose_in=False, sField_unique_id="cellid"
    )["area_min"]


def test_discrete_dominant_and_fractions(
    make_raster, make_grid_mesh, read_vector, tmp_path
):
    # 3 pixels class 1, 1 pixel class 2 -> mode 1, count 4, 75% / 25%
    values = np.array([[1, 1], [1, 2]], dtype=float)
    raster = make_raster(tmp_path / "d.tif", values=values, bbox=(0, 0, 1, 1))
    mesh = make_grid_mesh(
        tmp_path / "m.geojson", bbox=(0, 0, 1, 1), nrows=1, ncols=1
    )
    out = str(tmp_path / "out.geojson")
    extract.run_remap(
        out,
        mesh,
        [raster],
        _area_min(mesh),
        iFlag_remap_method_in=1,
        iFlag_discrete_in=True,
        iFlag_stat_in=False,
    )
    rows = read_vector(out)
    assert rows[1]["mode"] == 1
    assert rows[1]["count"] == 4
    assert rows[1]["percentage_1"] == pytest.approx(75.0, abs=1e-6)
    assert rows[1]["percentage_2"] == pytest.approx(25.0, abs=1e-6)
