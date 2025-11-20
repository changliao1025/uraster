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




affiliations:
 - name: Atmospheric, Climate, and Earth Sciences, Pacific Northwest National Laboratory, Richland, WA, USA
   index: 1


 - name: University of Calgary, Calgary, Canada
   index: 2

date: 13 Nov 2025


bibliography: paper.bib
---


# Summary


Converting existing **structured raster datasets** into an **unstructured mesh** is a necessary preliminary step for configuring numerical models that require spatially varying inputs and parameters. This repository offers the Python package, **`uraster`**, to facilitate this conversion. The package is highly flexible, designed to handle both continuous and categorical raster data, and incorporates options for **mass conservation** and various **interpolation methods** to ensure the resulting unstructured mesh accurately represents the original data. Furthermore, `uraster` integrates the industrial-standard 3D visualization library, **`GeoVista`**, allowing users to visualize the mesh data on a sphere, making it suitable for both regional and global scale applications. The package is intended to support variable-resolution unstructured mesh-based hydrologic and land surface models.


# Statement of need


There is an emergencing interest in using unstructured meshes for hydrologic and land surface modeling, as they can provide higher resolution in areas of interest while maintaining coarser resolution elsewhere, thus optimizing computational resources. Besides, unstructured meshes can better represent complex land surface features, such as rivers and lakes. However, most existing climate datasets are in the form of structured raster data, which cannot be directly used for unstructured mesh-based models. Therefore, there is a need for a tool that can convert structured raster datasets into unstructured meshes while preserving the integrity of the original data. The `uraster` package addresses this need by providing a flexible and efficient solution for this conversion process.


Through industry-standard geospatial libraries, `uraster` calculates the topological relationships (e.g., intersect and contain) between the raster datasets and unstructured meshes, performing various interpolation methods to extract raster information for each individual mesh cell. For example:


* For raster datasets with continuous values, when any GDAL-supported resampling method is configured, `uraster` extracts all pixels within or intersecting each mesh cell and performs geostatistical analysis.
* For raster datasets with continuous values, when the weighted area resampling method is configured (ensuring mass conservation), `uraster` calculates the area of each pixel within or overlapping the mesh cell and uses this information to compute a weighted average of the pixel values.
* For raster datasets with categorical values, `uraster` extracts all pixels that fall within or intersect each mesh cell. It then saves all identified categories along with their corresponding percentages to the mesh cell, effectively calculating the fractional representation of each category within that cell.


All the operations are performed using the standard GDAL/OGR APIs, ensuring robust data handling and compatibility with a wide range of geospatial data formats. The package also supports mesh cells that cross the International Date Line (IDL), making it suitable for global applications.


# Model features


- **GDAL-Native Vector Handling**: Uses the standard GDAL/OGR APIs for defining unstructured mesh cells, and performing projection-aware geospatial operations. It also supports mesh cells that cross the International Date Line (IDL).


- **Standard Vector I/O**: Instead of directly operating on various mesh standards, it utilizes standard geographic information system vector formats (e.g., GeoJSON) for mesh operations, ensuring broad compatibility. It supports major output formats through GDAL APIs.


- **Projection-Aware Operations**: Handles (raster datasets) map projection differences to ensure accurate aggregation of raster values within each mesh cell polygon.


- **Interactive GeoVista API**: Offers essential functions to visualize the input and output vector layers on a 3D sphere, including support for mesh cells that cross the International Date Line (IDL).


# State of the field


There are several existing tools that provide similar functionality, such as `raster2dggs`, `xESMF`, `pyresample`, and `rasterio`. However, these tools often have different design philosophies and may not fully meet the needs of hydrologic and land surface modeling users. For example,


* `raster2dggs` is designed for converting raster data into discrete global grid systems (DGGs), which does not support other types of unstructured meshes, such as the Model for Prediction Across Scales (MPAS) mesh and Triangulated irregular network (TIN) meshes.
* `xESMF` is primarily focused on regridding and interpolation of climate data that are quadrilateral and does not support non-quadrilateral meshes, such as the MPAS and TIN meshes.
* `pyresample` and `rasterio` are powerful tools for resampling and reprojecting raster data. However, neither package supports unstructured meshes. This limitation is partly due to pyresample's dependence on rasterio, which is primarily designed for handling structured raster data and does not have built-in functionality for unstructured mesh formats.


# Acknowledgment


The model described in this repository was supported by the following:


* the U.S. Department of Energy Office of Science Biological and Environmental Research through the Earth System Development program as part of the Energy Exascale Earth System Model (E3SM) project.


* the Earth System Model Development and Regional and Global Model Analysis program areas of the U.S. Department of Energy, Office of Science, Biological and Environmental Research program as part of the multi-program, collaborative Integrated Coastal Modeling (ICoM) project.


* the Earth System Model Development and Regional and Global Model Analysis program areas of the U.S. Department of Energy, Office of Science, Biological and Environmental Research program as part of the multi-program, collaborative Interdisciplinary Research for Arctic Coastal Environments (InteRFACE) project.


A portion of this research was performed using PNNL Research Computing at Pacific Northwest National Laboratory.


PNNL is operated for DOE by Battelle Memorial Institute under contract DE-AC05-76RL01830.


# References


Liao, C., & Li, M. (2025). URaster: A Python Package for Converting Structured Raster Data into Unstructured Meshes (v0.1.1). Zenodo. https://doi.org/10.5281/zenodo.17613497









