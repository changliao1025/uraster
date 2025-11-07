These examples demonstrate case studies that accept multiple raster inputs.

When using multiple rasters, ensure that the list of raster file paths is provided in the configuration dictionary under the key `'aFilename_source_raster'`.

In this scenario, all the raster files are considered during the GDAL WARP operation, even if they may overlap or have different spatial resolutions.

This is also useful if a global analysis is desired, where multiple raster datasets covering different regions are processed together, especially when Greenland and Antarctica DEMs are involved.

Because different raster datasets may have different spatial references (CRS), the GDAL WARP operation will reproject them to match the target mesh's CRS before performing zonal statistics.


