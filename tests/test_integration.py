import numpy as np
import os
import pytest
from uraster.classes.uraster import uraster


def test_end_to_end_continuous(make_raster, make_grid_mesh, read_vector, tmp_path):
    raster = make_raster(
        tmp_path / "e.tif", values=np.full((4, 4), 10.0), bbox=(0, 0, 1, 1)
    )
    mesh = make_grid_mesh(
        tmp_path / "m.geojson", bbox=(0, 0, 1, 1), nrows=1, ncols=1
    )
    out = str(tmp_path / "uraster.geojson")

    aConfig = {
        "sFilename_source_mesh": mesh,
        "aFilename_source_raster": [raster],
        "sFilename_target_mesh": out,
    }
    pRaster = uraster(aConfig)
    assert pRaster.setup() is True
    pRaster.run_remap()

    assert os.path.exists(out)
    rows = read_vector(out)
    assert rows[1]["mean"] == pytest.approx(10.0, abs=1e-6)
