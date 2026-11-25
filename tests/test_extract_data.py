import numpy as np
import pytest

from uraster import URaster, extract_raster_by_mesh, extract_raster_by_polygon


def test_extract_raster_by_mesh(make_raster, make_grid_mesh, tmp_path):
    values = np.arange(16, dtype=np.float32).reshape(4, 4)
    raster_path = make_raster(
        tmp_path / "test_extract.tif", values=values, bbox=(0, 0, 1, 1)
    )
    mesh_path = make_grid_mesh(
        tmp_path / "mesh_extract.geojson", bbox=(0, 0, 1, 1), nrows=2, ncols=2
    )

    # 1. Functional API
    data_dict = extract_raster_by_mesh(mesh_path, [raster_path])
    assert len(data_dict) == 4
    for cid, pixels in data_dict.items():
        assert isinstance(pixels, np.ndarray)
        assert len(pixels) > 0

    all_extracted = np.sort(np.concatenate(list(data_dict.values())))
    assert np.allclose(all_extracted, np.arange(16))

    # 2. Object-oriented URaster API
    pRaster = URaster(
        {
            "sFilename_source_mesh": mesh_path,
            "aFilename_source_raster": [raster_path],
            "sField_unique_id": "cellid",
        }
    )
    oo_data = pRaster.extract_raster_data()
    assert len(oo_data) == 4
    for cid in data_dict:
        assert np.array_equal(data_dict[cid], oo_data[cid])


def test_extract_raster_by_polygon(make_raster, tmp_path):
    values = np.full((4, 4), 42.0, dtype=np.float32)
    raster_path = make_raster(
        tmp_path / "test_poly.tif", values=values, bbox=(0, 0, 1, 1)
    )
    wkt = "POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))"

    pixels = extract_raster_by_polygon(wkt, [raster_path])
    assert isinstance(pixels, np.ndarray)
    assert len(pixels) == 16
    assert np.allclose(pixels, 42.0)
