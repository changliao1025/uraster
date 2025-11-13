##################
What is uraster?
##################

*********
Overview
*********

**uraster** is a Python package designed to convert or transfer structured raster datasets into unstructured mesh formats. It bridges the gap between structured raster data and unstructured mesh-based numerical models, leveraging GDAL/OGR for robust geospatial data handling.

The package addresses the common challenge in computational modeling where structured raster data needs to be mapped onto unstructured meshes for various scientific applications, including climate modeling, environmental simulation, and geophysical analysis.

uraster provides a seamless workflow to extract statistical information from raster datasets and assign these values to unstructured mesh cells through advanced zonal statistics operations.

***********
Development
***********

uraster is developed as an open-source project hosted on GitHub:
https://github.com/changliao1025/uraster

*********
Objective
*********

Many computational models use unstructured meshes to represent complex geometries and variable resolution requirements. However, input datasets are often available in structured raster formats. Traditional interpolation methods may not preserve important statistical properties or spatial relationships.

uraster addresses this gap by providing:

- **Projection-aware operations**: Handles coordinate system differences between raster and mesh data
- **Statistical preservation**: Maintains important statistical properties during the transfer process
- **Flexible mesh support**: Works with various unstructured mesh formats including MPAS, DGGS, and TIN
- **GDAL integration**: Leverages industry-standard geospatial libraries for robust data handling

*****************
Target Audience
*****************

uraster is designed for researchers, scientists, and developers working with:

- Climate and weather modeling
- Environmental simulation
- Geophysical analysis
- Computational fluid dynamics
- Any application requiring raster-to-mesh data conversion

Users should have basic familiarity with Geographic Information Systems (GIS) concepts, including raster and vector data, coordinate systems, and map projections.

*****************
Key Features
*****************

1. **GDAL-Native Vector Handling**: Uses standard GDAL/OGR for mesh operations
2. **Standard Vector I/O**: Supports common GIS formats (GeoJSON, Shapefile, etc.)
3. **Interactive Visualization**: Built-in 3D visualization capabilities with GeoVista
4. **Multiple Mesh Types**: Support for MPAS, DGGS, TIN, and custom unstructured meshes
