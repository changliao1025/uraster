from uraster import utility


def test_idl_crossing_polygon_is_split(make_polygon_mesh, read_vector, tmp_path):
    # ring spanning lon 170..190 (190 normalizes to -170) -> crosses the IDL
    ring = [(170, 0), (190, 0), (190, 10), (170, 10), (170, 0)]
    mesh = make_polygon_mesh(tmp_path / "idl.geojson", polygons=[ring])
    out = str(tmp_path / "idl_fixed.geojson")

    ok = utility.fix_mesh_longitude_range_and_idl_crossing(
        mesh, out, handle_idl_crossing=True
    )
    assert ok is True

    rows = read_vector(out)
    assert len(rows) == 1
    geom = rows[1]["__geom__"]
    assert geom.GetGeometryName() == "MULTIPOLYGON"
    assert geom.GetGeometryCount() == 2

    # every coordinate is within the valid longitude range
    for i in range(geom.GetGeometryCount()):
        poly = geom.GetGeometryRef(i)
        ring_geom = poly.GetGeometryRef(0)
        for p in range(ring_geom.GetPointCount()):
            x, _ = ring_geom.GetPoint_2D(p)
            assert -180.0 <= x <= 180.0


def test_non_crossing_polygon_preserved(make_polygon_mesh, read_vector, tmp_path):
    # a normal polygon must pass through as a single POLYGON (regression guard)
    ring = [(10, 0), (20, 0), (20, 10), (10, 10), (10, 0)]
    mesh = make_polygon_mesh(tmp_path / "plain.geojson", polygons=[ring])
    out = str(tmp_path / "plain_fixed.geojson")

    ok = utility.fix_mesh_longitude_range_and_idl_crossing(
        mesh, out, handle_idl_crossing=True
    )
    assert ok is True

    rows = read_vector(out)
    assert len(rows) == 1
    assert rows[1]["__geom__"].GetGeometryName() == "POLYGON"
