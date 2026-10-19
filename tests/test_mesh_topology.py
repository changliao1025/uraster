from uraster import utility


def test_rebuild_mesh_topology_counts(make_grid_mesh, tmp_path):
    path = make_grid_mesh(
        tmp_path / "grid.geojson", bbox=(0, 0, 3, 3), nrows=3, ncols=3
    )
    info = utility.rebuild_mesh_topology(
        path, iFlag_verbose_in=False, sField_unique_id="cellid"
    )
    assert info is not None
    assert info["num_cells"] == 9
    assert info["num_polygns"] == 9
    assert list(info["cell_ids"]) == list(range(1, 10))
    assert info["area_min"] > 0
