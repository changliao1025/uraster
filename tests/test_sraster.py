import numpy as np
import pytest
from uraster.classes.sraster import sraster


def test_metadata_matches_fixture(make_raster, tmp_path):
    path = make_raster(
        tmp_path / "meta.tif",
        values=np.zeros((4, 4)),
        bbox=(10.0, 20.0, 14.0, 24.0),  # 4x4 deg over 4x4 px -> 1 deg/px
        epsg=4326,
        nodata=-9999.0,
    )
    r = sraster(sFilename_in=path)
    r.read_metadata()
    assert r.ncolumn == 4
    assert r.nrow == 4
    assert r.dResolution_x == pytest.approx(1.0)
    assert r.dResolution_y == pytest.approx(1.0)
    assert r.dNoData == pytest.approx(-9999.0)
    assert r.pSpatialRef.GetAuthorityCode(None) == "4326"


def test_convert_to_wgs84_preserves_constant(make_raster, tmp_path):
    # constant field in Web Mercator; a small patch near the equator/prime meridian
    path = make_raster(
        tmp_path / "merc.tif",
        values=np.full((8, 8), 7.0),
        bbox=(0.0, 0.0, 10000.0, 10000.0),  # metres (EPSG:3857)
        epsg=3857,
        nodata=-9999.0,
    )
    r = sraster(sFilename_in=path)
    r.read_metadata()
    r2 = r.convert_to_wgs84()
    assert r2.pSpatialRef.GetAuthorityCode(None) == "4326"
    data = r2.read_data()
    valid = data[data != -9999.0]
    assert valid.size > 0
    assert np.allclose(valid, 7.0)
