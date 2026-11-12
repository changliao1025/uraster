import numpy as np
import pytest

from uraster.classes.sraster import SRaster


def test_metadata_matches_fixture(make_raster, tmp_path):
    path = make_raster(
        tmp_path / "meta.tif",
        values=np.zeros((4, 4)),
        bbox=(10.0, 20.0, 14.0, 24.0),  # 4x4 deg over 4x4 px -> 1 deg/px
        epsg=4326,
        nodata=-9999.0,
    )
    r = SRaster(sFilename_in=path)
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
    r = SRaster(sFilename_in=path)
    r.read_metadata()
    r2 = r.convert_to_wgs84()
    assert r2.pSpatialRef.GetAuthorityCode(None) == "4326"
    data = r2.read_data()
    valid = data[data != -9999.0]
    assert valid.size > 0
    assert np.allclose(valid, 7.0)


def test_resample_resolution_and_values(make_raster, tmp_path):
    path = make_raster(
        tmp_path / "src.tif",
        values=np.full((4, 4), 15.0),
        bbox=(0.0, 0.0, 4.0, 4.0),  # 1.0 deg/px
        epsg=4326,
        nodata=-9999.0,
    )
    r = SRaster(sFilename_in=path)
    r.read_metadata()
    assert r.dResolution_x == pytest.approx(1.0)

    # Resample to 0.25 deg/px (higher resolution)
    r_high = r.resample(dResolution_x=0.25, dResolution_y=0.25, sResampleAlg="near")
    assert r_high.dResolution_x == pytest.approx(0.25)
    assert r_high.dResolution_y == pytest.approx(0.25)
    assert r_high.ncolumn == 16
    assert r_high.nrow == 16
    data = r_high.read_data()
    assert np.allclose(data, 15.0)
