# Example data

To keep the code repository lightweight, `uraster` separates code from data. This
`data/` directory holds only the **outputs** produced by the example scripts; the
**inputs** (meshes and rasters) are distributed separately.

## Directory layout

```
data/
  example_N/
    input/      # source mesh + raster(s) — populated on demand (see below)
    output/     # results committed to the repository
```

Each `example_N/output/` folder contains the products of the corresponding example:

- `uraster.geojson` — the target mesh with the remapped variable.
- `mesh.jpg` — 3D visualization of the source mesh on a sphere.
- `uraster.png` — 3D visualization of the remapped variable.
- `uraster.mp4` / `.gif` — optional rotation animation.

## Input data

Because spatial datasets exceed GitHub's file-size limits, input meshes and rasters are
hosted in a separate repository and served through its GitHub releases:

**https://github.com/changliao1025/uraster_data**

The example scripts download and cache these inputs automatically on first run using
[Pooch](https://www.fatiando.org/pooch/), so you normally do not need to fetch them by
hand. See `get_example_paths()` in `uraster/utility.py` and the caching helpers
(`get_cache_location()`, `clear_cache()`) for details.

See the README in the `uraster_data` repository for a full description of each input
dataset. A catalog of all examples — mesh type, raster, data type, and remapping method
— is maintained in [`examples/readme.md`](../examples/readme.md).
