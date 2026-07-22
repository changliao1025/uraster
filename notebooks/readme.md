# URaster Notebooks

This directory contains Jupyter notebook versions of the scripts in
[`examples`](../examples). Each notebook reproduces one example end to end —
downloading inputs, running the remap, and rendering the 3D visualizations inline — so
you can explore `uraster` interactively without leaving the browser.

## Running the notebooks

You can launch the notebooks in the cloud with **Binder** (no local install required):

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/changliao1025/uraster/main?urlpath=%2Fdoc%2Ftree%2Fnotebooks%2Fexample_1%2Frun_rhpx_example.ipynb)

To run locally, install `uraster` (see the top-level [README](../README.md)) along with
the notebook extras, then open any notebook with Jupyter:

```bash
pip install "uraster[notebook]"
jupyter lab notebooks/example_1/run_rhpx_example.ipynb
```

3D rendering uses [GeoVista](https://geovista.readthedocs.io/) (PyVista/VTK). In headless
environments the notebooks detect the absence of a display and fall back to off-screen
rendering automatically.

## Notebook catalog

Each notebook mirrors the example of the same number. See
[`examples/readme.md`](../examples/readme.md) for the full description of each case
(mesh type, raster data, data type, and remapping method).

| # | Notebook | Mesh type | Highlights |
|---|----------|-----------|-----------|
| 1 | `example_1/run_rhpx_example.ipynb` | rHEALPix (global) | continuous raster, average |
| 2 | `example_2/run_isea3h_example.ipynb` | ISEA3H | categorical raster, dominant class |
| 3 | `example_3/run_rhpx_example.ipynb` | rHEALPix (China) | continuous raster, average |
| 4 | `example_4/run_rhpx_example.ipynb` | rHEALPix (China, finer) | finer mesh resolution |
| 5 | `example_5/run_isea7h_example.ipynb` | ISEA7H | multi-tile raster mosaic |
| 6 | `example_6/run_isea7h_example.ipynb` | ISEA7H (high latitude) | multi-tile, mixed CRS, polar |
| 7 | `example_7/run_isea3h_example.ipynb` | ISEA3H | multi-tile, mixed CRS, categorical |
| 8 | `example_8/run_tin_example.ipynb` | TIN | triangulated irregular network mesh |
| 9 | `example_9/run_mpas_example.ipynb` | MPAS | global MPAS mesh |
| 10 | `example_10/run_rhpx_example.ipynb` | rHEALPix (China) | mass-conserving weighted average |
