# URaster Examples

This directory contains runnable Python scripts that demonstrate `uraster` across a
range of unstructured mesh types, raster data types, and remapping methods. Each
example converts a structured raster dataset onto an unstructured mesh and produces a
target mesh (with the remapped variable) plus static and animated 3D visualizations.

Jupyter notebook versions of these scripts are available in the [`notebooks`](../notebooks)
directory.

## How the examples work

Every example follows the same workflow:

1. **Download inputs** with [Pooch](https://www.fatiando.org/pooch/) via
   `get_example_paths(example_number=N)`. Input meshes and rasters are hosted in the
   companion data repository, [`uraster_data`](https://github.com/changliao1025/uraster_data),
   and cached locally on first run.
2. **Configure** the run through an `aConfig` dictionary:
   - `sFilename_source_mesh` — the input unstructured mesh (vector).
   - `aFilename_source_raster` — a list of one or more source raster tiles.
   - `sFilename_target_mesh` — the output mesh path.
   - `iFlag_discrete` — set to `1` for categorical rasters (optional).
3. **Run** the pipeline:
   ```python
   pRaster = uraster(aConfig)
   pRaster.setup()
   pRaster.report_inputs()
   pRaster.visualize_source_mesh(sFilename_out="mesh.jpg", ...)
   pRaster.run_remap(...)
   pRaster.visualize_target_mesh(sFilename_out="uraster.png", ...)
   pRaster.cleanup()
   ```

Outputs are written to `data/example_N/output/` relative to the working directory
(`mesh.jpg`, `uraster.geojson`, `uraster.png`, and an `.mp4`/`.gif` animation).

### Remapping methods

`run_remap()` selects the aggregation strategy through its arguments:

| Method | How to invoke | Used in |
|--------|---------------|---------|
| Average (default) | `run_remap()` | 1, 3, 4, 5, 6, 8, 9 |
| Dominant (categorical) | `aConfig["iFlag_discrete"] = 1` and `run_remap(iFlag_discrete_in=True)` | 2, 7 |
| Weighted (mass-conserving) | `run_remap(iFlag_weighted_average_in=True)` | 10 |

## Example catalog

| # | Script | Mesh type | Mesh file | Cell size | Raster data | Tiles | Data type | Method |
|---|--------|-----------|-----------|-----------|-------------|-------|-----------|--------|
| 1 | `example_1/run_rhpx_example.py` | rHEALPix (global) | `rhealpix_global_res3.geojson` | 116,744 km² | EDGAR CH₄ manure-management inventory, 2015 (ton) | single | continuous | Average |
| 2 | `example_2/run_isea3h_example.py` | ISEA3H (Toronto & Buffalo) | `isea3h_bbox_res15.geojson` | 3.55 km² | World Settlement Footprint 2015 (255/0 binary) | single | discrete | Dominant |
| 3 | `example_3/run_rhpx_example.py` | rHEALPix (mainland China) | `rhealpix_China_res3.geojson` | 116,744 km² | Coal-mine methane emissions, China 2020 (Mg km⁻² a⁻¹) | single | continuous | Average |
| 4 | `example_4/run_rhpx_example.py` | rHEALPix (mainland China, finer) | `rhealpix_China_res6.geojson` | 160 km² | same as example 3 | single | continuous | Average |
| 5 | `example_5/run_isea7h_example.py` | ISEA7H (inland eastern US) | `isea7h_bbox_res5.geojson` | 3,034.84 km² | Hansen Global Forest Change tree cover, 2000 | multi-tile (4) | continuous | Average |
| 6 | `example_6/run_isea7h_example.py` | ISEA7H (Nunavut, high latitude) | `isea7h_bbox_res8.geojson` | 1.26 km² | ArcticDEM + FABDEM elevation | multi-tile (9, mixed CRS) | continuous | Average |
| 7 | `example_7/run_isea3h_example.py` | ISEA3H (Rocky Mountains) | `isea3h_bbox_res12.geojson` | 95.98 km² | Esri land use/land cover from Sentinel-2, 2024 | multi-tile (4, mixed CRS) | discrete | Dominant |
| 8 | `example_8/run_tin_example.py` | TIN (Susquehanna basin, JIGSAW) | `tin.geojson` | — | HydroSHEDS global DEM | single | continuous | Average |
| 9 | `example_9/run_mpas_example.py` | MPAS (global) | `mpas.geojson` | — | same as example 8 | single | continuous | Average |
| 10 | `example_10/run_rhpx_example.py` | rHEALPix (mainland China) | same as example 3 | — | same as example 3 (resampled first) | single | continuous | Weighted average |

## What the set demonstrates

- **Mesh variety** — rHEALPix, ISEA3H, and ISEA7H discrete global grids, plus TIN and
  MPAS meshes, showing `uraster` is not limited to DGGS.
- **Continuous vs. categorical rasters** — examples 2 and 7 use the categorical
  (dominant/fractional) path; the rest are continuous.
- **Single vs. multi-tile inputs** — examples 5, 6, and 7 mosaic multiple raster tiles;
  examples 6 and 7 additionally mix coordinate reference systems, exercising on-the-fly
  reprojection to WGS84.
- **Aggregation methods** — nearest/average vs. mass-conserving weighted average
  (example 10), including a resample-then-remap workflow.
- **Spatial scale** — global (examples 1, 9), regional (examples 3–7), and
  high-latitude/polar (example 6).

## Data sources

| # | Source | Reference |
|---|--------|-----------|
| 1 | [EDGAR](https://data.jrc.ec.europa.eu/dataset/b54d8149-2864-4fb9-96b9-5fd3a020c224) | Crippa, M., et al. (2021). High resolution temporal profiles in the Emissions Database for Global Atmospheric Research (EDGAR). |
| 2 | [WSF2015](https://download.geoservice.dlr.de/WSF2015/) | Marconcini, M., et al. (2020). Outlining where humans live, the World Settlement Footprint 2015. *Scientific Data*, 7(1), 242. |
| 3, 4, 10 | [Coal mine methane, China](https://pubs.acs.org/doi/full/10.1021/acs.estlett.9b00294) | Sheng, J., et al. (2019). Bottom-Up Estimates of Coal Mine Methane Emissions in China. *Environmental Science & Technology Letters*, 6(8), 473–478. |
| 5 | [Hansen Global Forest Change](https://earthenginepartners.appspot.com/science-2013-global-forest/download_v1.7.html) | Hansen, M. C., et al. (2013). High-Resolution Global Maps of 21st-Century Forest Cover Change. *Science*, 342(6160), 850–853. |
| 6 | [FABDEM](https://data.bris.ac.uk/data/dataset/25wfy0f9ukoge2gs7a5mqpq2j7), [ArcticDEM](https://www.arcgis.com/apps/webappviewer/index.html?id=aff5fa8f5d5548c6bff44cc8be385f61) | Hawker, L., et al. (2022). A 30 m global map of elevation with forests and buildings removed. *Environmental Research Letters*, 17(2), 024016. Porter, C., et al. (2023). ArcticDEM — Mosaics, Version 4.1, Harvard Dataverse. doi:10.7910/DVN/3VDC4W |
| 7 | [Esri Land Cover Explorer](https://livingatlas.arcgis.com/landcoverexplorer/) | Karra, K., et al. (2021). Global land use/land cover with Sentinel-2 and deep learning. *IGARSS 2021*. |
| 8, 9 | [HydroSHEDS](https://www.hydrosheds.org/) | Lehner, B., et al. (2008). New global hydrography derived from spaceborne elevation data. *Eos, Transactions AGU*, 89(10), 93–94. |

## Running an example

```bash
cd examples/example_1
python run_rhpx_example.py
```

On first run the required inputs are downloaded to the Pooch cache; subsequent runs
reuse the cached files.
