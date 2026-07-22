# URaster Offline Test Suite — Design

**Date:** 2026-07-21
**Status:** Approved (design)
**Context:** JOSS submission (issue openjournals/joss-reviews#10556) was desk-rejected.
Editor feedback: unit tests are too minimal, and it is not obvious how a user would
confirm the software works correctly (no reference results to compare against; docs not
compiled). This spec covers the **tests** portion only. Documentation/ReadTheDocs is a
separate later phase.

## Goal

Provide an automated, **offline**, correctness-checking test suite that:

1. Runs without any network access, on Linux/Windows/macOS, under `pytest`.
2. Verifies the numerical claims made in `paper.md` — continuous averaging, categorical
   composition, mass-conserving weighted aggregation, multi-tile mosaicking, on-the-fly
   reprojection, and antimeridian (IDL) handling — by comparing against **hand-computable
   expected values** (inline analytic assertions, not opaque golden files).
3. Is wired into CI so a reviewer sees a green badge, and documents how to run it.

Success criterion: `pytest` (default markers) passes offline in the `uraster` conda env
and in CI, and each core code path asserts a known-correct numeric result.

## Non-goals (this round)

- Visualization / GeoVista rendering correctness. Only an import-level smoke check;
  rendering is heavy and non-deterministic.
- Golden reference files. We deliberately chose inline analytic assertions.
- Documentation / ReadTheDocs build. Separate phase.
- 100% line coverage. We target the scientifically meaningful paths, not every branch.

## Environment

- The package requires `numpy`, `gdal>=3.0.0`, `pooch`, `pyearth>=0.2.3`,
  `pyearthviz3d>=0.1.2`; dev tools add `pytest`, `pytest-cov`.
- No ready-to-use env currently exists on the dev machine. The `uraster` env is created
  from `environment.yml`, then `pytest`/`pytest-cov` added and the package installed with
  `pip install -e .`. Environment/package installation is performed by the user, not by
  the agent.
- Tests target the `uraster` conda env.

## Test framework and level

- **Framework:** `pytest` with fixtures (`tmp_path`, `tmp_path_factory`, `approx`,
  `parametrize`, markers). Existing `unittest`-style files continue to run under pytest.
- **Level:** both — unit tests against the low-level functions
  (`operation.extract.run_remap`, `operation.intersect.run_remap`,
  `utility.rebuild_mesh_topology`, `sraster` methods) for deterministic math, **plus** one
  integration test through the public `uraster` class (`setup()` → `run_remap()`, with
  visualization skipped) to prove the real user path.

## Relevant implementation facts (verified against source)

- `uraster.run_remap` dispatch (`uraster/classes/uraster.py`):
  - `iFlag_weighted_average_in=True` → `intersect.run_remap` (builds a pixel-footprint
    "raster mesh" via `sraster.create_raster_mesh()`, then area-weights).
  - otherwise → `extract.run_remap(..., iFlag_remap_method_in in {1,2,3})`.
  - `iFlag_discrete_in=True` forces method 1 and disables extra statistics.
- `extract.run_remap(sFilename_target_mesh, sFilename_source_mesh,
  aFilename_source_raster, dArea_min, iFlag_remap_method_in=1, iFlag_stat_in=True,
  iFlag_discrete_in=False, ..., sField_unique_id="cellid")`.
- Output vector schema (`extract.py`):
  - Continuous: fields `cellid`, `area`, `mean` (+ `min`, `max`, `std` when
    `iFlag_stat_in=True`).
  - Discrete: fields `cellid`, `area`, `mode` (dominant class as int, `-1` if none),
    `count`, and one `percentage_{int(value)}` field per unique raster category.
- `utility.rebuild_mesh_topology(sFilename_mesh, iFlag_verbose_in, sField_unique_id)`
  returns a dict including `num_cells`, `num_polygns`, `cell_ids` (array), and area stats
  (`area_min`, ...). `dArea_min` from this feeds `extract.run_remap`.

## File layout

```
tests/
  conftest.py               # shared fixtures + GDAL fixture-builders
  test_basic.py             # KEEP as-is (import/instantiation smoke)
  test_sraster.py           # metadata, CRS, geotransform, reprojection
  test_mesh_topology.py     # rebuild_mesh_topology determinism
  test_remap_continuous.py  # average + min/max/std
  test_remap_weighted.py    # mass-conserving area weighting (intersect)
  test_remap_discrete.py    # dominant class + percentage fractions
  test_remap_multitile.py   # tile mosaic + on-the-fly reprojection
  test_remap_idl.py         # antimeridian-crossing cell
  test_integration.py       # end-to-end via uraster class (no visualization)
  test_workflow.py          # EXISTING network test — marked, bug-fixed, skipped by default
```

Rationale for the split: one file per concern keeps each file small and focused, matching
the reviewer's readability concern. Files may be merged during implementation only if a
concern turns out trivial; the concern coverage below is the contract, not the file count.

## Fixture builders (`conftest.py`)

All fixtures write to `tmp_path` — nothing is committed to the repo. Helpers:

- `make_raster(path, values, bbox, crs="EPSG:4326", nodata=None, dtype)` — writes a small
  GeoTIFF whose pixel values are supplied as a 2-D array (or constant), with a known
  geotransform derived from `bbox` and array shape.
- `make_grid_mesh(path, bbox, nrows, ncols, crs="EPSG:4326", id_field="cellid")` — writes a
  GeoJSON polygon mesh of `nrows×ncols` rectangular cells with sequential `cellid`
  starting at 1 (reusing the pattern already in the existing `create_test_mesh`).
- `make_polygon_mesh(path, rings, crs, id_field="cellid")` — writes an explicit-geometry
  mesh for special cases (single cell, IDL-crossing cell).

Values are always chosen so the expected aggregate is computable by hand.

## Correctness anchors

Each row is an independent test with an inline expected value. Tolerances: `abs=1e-6` for
exact-constant cases; a documented looser tolerance (~`1e-3` relative) where sub-pixel
resampling or reprojection introduces small error.

| # | Test file | Fixture | Assertion |
|---|-----------|---------|-----------|
| 1 | test_sraster | 4×4 raster, known geotransform + EPSG:4326 + nodata | `ncolumn==4`, `nrow==4`, resolution, `pSpatialRef` EPSG 4326, `dNoData` match |
| 2 | test_sraster | raster in EPSG:32617 | `convert_to_wgs84()` produces a WGS84 raster; extent transforms sanely |
| 3 | test_mesh_topology | 3×3 grid mesh | `num_cells==9`, `num_polygns==9`, `cell_ids==1..9`, `area_min>0` |
| 4 | test_remap_continuous | raster all = 10.0 over single-cell mesh | `mean == approx(10.0)` |
| 5 | test_remap_continuous | raster half 0 / half 20 over single cell | `mean == approx(10.0)` |
| 6 | test_remap_continuous | multi-cell grid, each cell over a constant region | each cell `mean` == its region constant |
| 7 | test_remap_continuous | values 0..15, `iFlag_stat_in=True` | `min`/`max`/`std` match numpy on the same array |
| 8 | test_remap_weighted | known pixel values + cell of known fractional pixel overlap | area-weighted mean == hand-computed value; total mass (Σ value·area) conserved within tol |
| 9 | test_remap_discrete | cell = 3 px class 1 + 1 px class 2 | `mode==1`, `count==4`, `percentage_1==approx(75.0)`, `percentage_2==approx(25.0)` |
| 10 | test_remap_multitile | two adjacent tiles, same constant value, one covering cell | `mean == approx(constant)` (mosaic works) |
| 11 | test_remap_multitile | projected-CRS raster (constant value) + WGS84 mesh | `mean == approx(constant)` (reprojection works) |
| 12 | test_remap_idl | single cell spanning ±180° over a raster with known values | correct `mean`, no crash, cell present in output |
| 13 | test_integration | raster + grid mesh via `uraster()` class | `setup()` then `run_remap()` (continuous) writes output GeoJSON whose `mean` matches expected; visualization not called |

Notes:
- #8 (mass conservation) and #12 (antimeridian) are the two headline paper claims and the
  fiddliest fixtures; they get explicit, commented derivations of the expected number.
- #12 stays a unit test against `extract.run_remap` (and/or the
  `utility.fix_mesh_longitude_range_and_idl_crossing` helper); it is not required to also
  be an integration test.

## Network test handling

- Register a `network` marker (in `pyproject.toml` `[tool.pytest.ini_options]`) and set
  `addopts = -m "not network"` so the default and CI runs skip downloads.
- Mark `tests/test_workflow.py` classes/functions with `@pytest.mark.network`.
- Fix its `setUpClass` bug: it currently builds fixtures with `iFlag_download_in=False`, so
  the raster files are never created and the remap tests fail. When run under `-m network`
  it must actually download; otherwise it is skipped. Also add real assertions there over
  time (not required for this round beyond making it non-broken and skipped).
- Move its heavy imports (`pystac_client`, `planetary_computer`, `elevation`) inside the
  functions or guard them, so collecting the default suite does not require those packages.

## CI

- `.github/workflows/ci.yml` currently triggers on `main`/`develop`; add `development`
  (the active branch) to `push`/`pull_request`.
- The default `pytest tests/ -v --cov=uraster` run already excludes the `network` marker
  via `addopts`; confirm the network-only imports are not collected.
- Keep the coverage upload.

## Risks / open questions (to resolve during implementation)

- The `intersect.run_remap` (weighted) output field name(s) for the aggregated value must
  be confirmed from source when writing test #8 (continuous `mean` vs a weighted field).
- Antimeridian fixture: confirm whether `extract.run_remap` expects the crossing cell in
  raw ±180 coordinates or normalized; build the fixture to match the documented behavior
  and, if needed, exercise `fix_mesh_longitude_range_and_idl_crossing` directly.
- Parallel path: `extract.run_remap` switches to multiprocessing above
  `iFeature_parallel_threshold` (10000 features). All fixtures stay tiny (well under the
  threshold) so tests run single-process and deterministically.

## Deliverables

1. `tests/conftest.py` with fixture builders.
2. The test files listed above with the anchors implemented.
3. `pyproject.toml` pytest config (`markers`, `addopts`).
4. `ci.yml` trigger fix.
5. A short "Running the tests" section for the README (wording finalized in the docs phase,
   but the command is `pytest` in the `uraster` env).
