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


affiliations:
 - name: Atmospheric, Climate, and Earth Sciences, Pacific Northwest National Laboratory, Richland, WA, USA
   index: 1
date: 13 Nov 2025

bibliography: paper.bib
---

# Summary

Converting existing **structured raster datasets** into an **unstructured mesh** is a necessary preliminary step for configuring numerical models that require spatially varying inputs and parameters. This repository offers the Python package, **`uraster`**, to facilitate this conversion. The package is highly flexible, designed to handle both continuous and categorical raster data, and incorporates options for **mass conservation** and various **interpolation methods** to ensure the resulting unstructured mesh accurately represents the original data. Furthermore, `uraster` integrates the industrial-standard 3D visualization library, **`GeoVista`**, allowing users to visualize the mesh data on a sphere, making it suitable for both regional and global scale applications. The package is intended to support variable-resolution unstructured mesh-based hydrologic and land surface models.

# Statement of need

There is an emergencing interest in using unstructured meshes for hydrologic and land surface modeling, as they can provide higher resolution in areas of interest while maintaining coarser resolution elsewhere, thus optimizing computational resources. Besides, unstructured meshes can better represent complex land surfaace features, such as river and lake. However, most existing climate dataset are in the form of structured raster data, which cannot be directly used for unstructured mesh-based models. Therefore, there is a need for a tool that can convert structured raster datasets into unstructured meshes while preserving the integrity of the original data. The `uraster` package addresses this need by providing a flexible and efficient solution for this conversion process.


# Model features



# State of the field



# Acknowledgment

The model described in this repository was supported by the following:

* the U.S. Department of Energy Office of Science Biological and Environmental Research through the Earth System Development program as part of the Energy Exascale Earth System Model (E3SM) project.

* the Earth System Model Development and Regional and Global Model Analysis program areas of the U.S. Department of Energy, Office of Science, Biological and Environmental Research program as part of the multi-program, collaborative Integrated Coastal Modeling (ICoM) project.

* the Earth System Model Development and Regional and Global Model Analysis program areas of the U.S. Department of Energy, Office of Science, Biological and Environmental Research program as part of the multi-program, collaborative Interdisciplinary Research for Arctic Coastal Environments (InteRFACE) project.

A portion of this research was performed using PNNL Research Computing at Pacific Northwest National Laboratory.

PNNL is operated for DOE by Battelle Memorial Institute under contract DE-AC05-76RL01830.

# References

