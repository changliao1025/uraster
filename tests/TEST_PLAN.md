# URaster Test Improvement Plan

## Current status
- Existing tests in `tests/test_basic.py` only cover basic imports and simple instantiation of `uraster` and `sraster`.
- No regression or functional coverage for raster metadata, mesh topology, or remap workflows.
- The `uraster` conda environment exists at `C:\Users\chang_rnrigc1\miniconda3\envs\uraster`, but `pytest` is not installed there.
- The package itself is importable in the environment, and core helper utilities exist in `uraster/utility.py`.

## Findings
- `sraster` has runnable methods for `read_metadata()`, `read_data()`, and `create_raster_mesh()`.
- `uraster` supports workflow steps: `setup()`, `check_raster_files()`, `check_mesh_file()`, `run_remap()`, `rebuild_mesh_topology()`, and `cleanup()`.
- `utility.rebuild_mesh_topology()` is a strong target for deterministic testing with synthetic vector input.
- The existing examples use `pooch` to download data, but unit tests should avoid external network dependencies.

## Suggested next steps
1. Add a deterministic unittest that creates a tiny in-memory or temporary GeoJSON mesh and a small GDAL raster file.
2. Use `sraster` to read metadata and verify raster shape / CRS information.
3. Use `utility.rebuild_mesh_topology()` or `uraster.check_mesh_file()` to validate mesh topology output.
4. Add a regression test for `uraster.run_remap()` if a minimal synthetic raster/mesh combination can be created reliably.
5. Install `pytest` into the conda env or continue using `unittest` for test execution.

## Notes for later
- When resuming, verify the conda environment again and install missing dev packages if needed.
- Prefer tests that compare numerical outputs or metadata against fixed expected values rather than visual outputs.
- Consider adding a `tests/data/` subfolder with small generated test files for reproducible workflows.
