---
title: 'uraster: Structured Raster to Unstructured Mesh'


tags:
 - Python
 - raster
 - mesh
 - geographic information system


authors:
 - name: Chang Liao
   orcid: 0000-0002-7348-8858
   affiliation: 1


 - name: Mingke Li
   orcid: 0000-0001-6310-4964
   affiliation: 2

 - name: Bill Little
   orcid: 0000-0002-1345-9465
   affiliation: 3




affiliations:
 - name: Atmospheric, Climate, and Earth Sciences, Pacific Northwest National Laboratory, Richland, WA, USA
   index: 1


 - name: University of Calgary, Calgary, Canada
   index: 2

 - name: Met Office, Exeter, UK
   index: 3

date: 13 Nov 2025


bibliography: paper.bib
---


# Summary


Converting **structured raster datasets** into **unstructured meshes** is a common prerequisite for configuring numerical models that require spatially varying inputs and parameters. This repository offers the Python package, **`uraster`**, to facilitate this conversion. The package is highly flexible, designed to support both continuous and categorical raster data, and provides configurable options for **mass conservation** and multiple **interpolation methods** to ensure the resulting unstructured mesh accurately represents the original data. Furthermore, `uraster` integrates the industrial-standard 3D visualization library, **`GeoVista`**, enabling intuitive and interactive 3D visualization of mesh-based data on a sphere. These capabilities make `uraster` well-suited for regional to global scale applications and for variable resolution unstructured mesh-based hydrologic and land surface modeling.


# Statement of need


There is an emerging interest in using unstructured meshes for hydrologic and land surface modeling, as they provide higher resolution in areas of interest while maintaining coarser resolution elsewhere, thus optimizing computational resources. Unstructured meshes also provide greater flexibility for representing complex land surface features, such as rivers, lakes, and coastlines. However, most existing climate and environmental datasets are distributed as structured raster data, which cannot be directly ingested by unstructured mesh-based models. Therefore, it creates a clear need for a robust and general tool that can convert structured raster datasets into unstructured meshes while preserving the integrity and physical meaning of the original data.


The `uraster` package addresses this need by providing a flexible and efficient solution for raster-to-mesh conversion. Through industry-standard geospatial libraries, `uraster` calculates the topological relationships (e.g., intersect and contain) between the raster datasets and unstructured meshes, performing various interpolation methods to extract raster information for each mesh cell. For example:


* For raster datasets with continuous values, when a GDAL-supported resampling method is configured, `uraster` extracts all pixels within or intersecting each mesh cell and performs geostatistical analysis.
* For raster datasets with continuous values, when the weighted area resampling method is configured (ensuring mass conservation), `uraster` calculates the area of each pixel that overlaps a mesh cell and uses these areas as weights to compute aggregated values.
* For raster datasets with categorical values, `uraster` extracts all pixels that fall within or intersect each mesh cell and records the fractional contribution of each category, resulting in a quantitative representation of category composition at the mesh cell level.


All the operations are implemented using the standard GDAL/OGR APIs, ensuring robust data handling and compatibility with a wide range of geospatial data formats. The package also supports mesh cells that cross the antimeridian, enabling seamless application to regional and global scale modeling workflows.



# Model features

- **GDAL-Native Vector Handling**: Uses the standard GDAL/OGR APIs to define unstructured mesh cells and to perform projection-aware geospatial operations. The model supports mesh cells that cross the antimeridian, enabling robust global scale workflows.


- **Standard Vector I/O**: Utilizes standard geographic information system vector formats (e.g., GeoJSON) for mesh operations instead of operating on various model-specific mesh standards. Input and output operations are handled through GDAL APIs, ensuring broad compatibility with widely used geospatial data formats.


- **Projection-Aware Operations**: Handles differences in map projections between raster datasets and mesh geometries, ensuring accurate aggregation of raster values within each mesh cell polygon.


- **Interactive GeoVista API**: Offers essential functions to visualize the input and output vector layers on a 3D sphere, including support for mesh cells that cross the antimeridian and projecting an unstructured mesh to a planar coordinate reference system.


# State of the field


There are several existing tools that provide similar functionality, such as `raster2dggs`, `xESMF`, `pyresample`, and `rasterio`. While these packages are widely used and well maintained, they are built around different design assumptions and do not fully address the requirements of hydrologic and land surface modeling workflows that rely on general unstructured meshes. For example:


* `raster2dggs` is designed specifically for converting raster data into Discrete Global Grid Systems (DGGS), which does not support other types of unstructured meshes, such as the Model for Prediction Across Scales (MPAS) mesh and Triangulated Irregular Network (TIN) meshes.
* `xESMF` is primarily focused on regridding and interpolation of climate data that are quadrilateral and does not support non-quadrilateral meshes, such as the MPAS and TIN meshes.
*  `pyresample` and `rasterio` provide powerful capabilities for raster resampling, reprojection, and data access. However, neither package supports unstructured mesh representations. This limitation is partly due to pyresample's dependence on rasterio, which is primarily designed for handling structured raster data and does not have built-in functionality for unstructured mesh formats.


In contrast to these tools, `uraster`  is designed explicitly to bridge structured raster datasets and general unstructured meshes. By operating directly on vector-based mesh geometries and supporting both continuous and categorical data with optional mass-conserving aggregation, `uraster`  addresses a gap in existing geospatial software ecosystems for unstructured mesh-based environmental modeling.


# Acknowledgment


The model described in this repository was supported by the following:


* the U.S. Department of Energy Office of Science Biological and Environmental Research through the Earth System Development program as part of the Energy Exascale Earth System Model (E3SM) project.


* the Earth System Model Development and Regional and Global Model Analysis program areas of the U.S. Department of Energy, Office of Science, Biological and Environmental Research program as part of the multi-program, collaborative Integrated Coastal Modeling (ICoM) project.


* the Earth System Model Development and Regional and Global Model Analysis program areas of the U.S. Department of Energy, Office of Science, Biological and Environmental Research program as part of the multi-program, collaborative Interdisciplinary Research for Arctic Coastal Environments (InteRFACE) project.


A portion of this research was performed using PNNL Research Computing at Pacific Northwest National Laboratory.


PNNL is operated for DOE by Battelle Memorial Institute under contract DE-AC05-76RL01830.


# References


Liao, C., & Li, M. (2025). URaster: A Python Package for Converting Structured Raster Data into Unstructured Meshes (v0.1.1). Zenodo. https://doi.org/10.5281/zenodo.17613497









