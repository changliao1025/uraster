###########################
Frequently Asked Questions
###########################

**Q: What mesh formats does uraster support?**

A: uraster supports various unstructured mesh formats including:

- GeoJSON
- Shapefile
- MPAS NetCDF meshes
- DGGS (Discrete Global Grid System) meshes
- TIN (Triangulated Irregular Network) meshes

**Q: What raster formats are supported?**

A: uraster leverages GDAL, so it supports all GDAL-compatible raster formats including:

- GeoTIFF
- NetCDF
- HDF5
- ESRI ASCII Grid
- And many others

**Q: How does uraster handle coordinate system differences?**

A: uraster automatically reprojects data as needed using GDAL's projection capabilities. It ensures that raster and mesh data are in compatible coordinate systems before processing.

**Q: Can I use custom statistical functions?**

A: Currently, uraster provides standard zonal statistics (mean, sum, min, max, etc.). Custom functions may be supported in future versions.

**Q: Is uraster suitable for large datasets?**

A: Yes, uraster is designed to handle large-scale datasets efficiently. It uses chunking and memory management techniques for processing large raster files.

**Q: How do I report bugs or request features?**

A: Please use the GitHub issue tracker at https://github.com/changliao1025/uraster/issues