[![DOI](https://zenodo.org/badge/1060295281.svg)](https://doi.org/10.5281/zenodo.17613497)
[![PyPI version](https://badge.fury.io/py/uraster.svg)](https://badge.fury.io/py/uraster)
[![License](https://img.shields.io/badge/license-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![geovista](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/bjlittle/geovista/main/docs/assets/badge/v0.json)](https://geovista.readthedocs.io/)

# URaster: Structured Raster to Unstructured Mesh

## Overview

**URaster** is a Python package to convert structured raster datasets into unstructured mesh structure, designed to bridge the gap between structured raster datasets and unstructured mesh-based hydrologic and land surface models. It leverages GDAL/OGR for robust data handling.

## ✨ Core Features

- **GDAL-Native Vector Handling**: Uses the standard GDAL/OGR engine for defining unstructured mesh cells, and performing projection-aware geospatial operations. It also support mesh cells that cross the International Date Line (IDL).

- **Standard Vector I/O**: Instead of directly operating on various mesh standards, it utilizes standard geographic information system vector formats (e.g., GeoJSON) for mesh operations, ensuring broad compatibility. It supports transformation APIs between existing meshes and standard vector formats.

- **Projection-Aware Operations**: Handles (raster dateaset) map projection differences to ensure accurate aggregation of raster values within each polygon.

- **Interactive GeoVista API**: Offers simple functions to visualize the input and the output vector layers on a 3D sphere.

Static visualization of the unstructured mesh on a sphere using GeoVista:
![Unstructured Mesh Visualization](docs/figures/mesh.jpg)

Animated visualization of the unstructured raster on a sphere using GeoVista:
![Unstructured Raster Visualization](docs/figures/global_uraster.gif)

## 💻 Installation

URaster requires GDAL for vector handling and GeoVista (which relies on PyVista/VTK) for 3D visualization.

> ⚠️ **GDAL Note**: Installing GDAL's Python bindings can be complex via pip due to platform dependencies. We strongly recommend using Conda for a stable installation of GDAL and all dependencies.

### Install via Conda (Recommended)
```bash
# Create a new conda environment (recommended)
conda create -n uraster-env python=3.9
conda activate uraster-env

# Install uraster and all dependencies via conda
conda install -c conda-forge uraster
```

## 🚀 Quick Start

[Quickstart documentation](https://uraster.readthedocs.io/en/latest/quickstart.html)

Example datasets are provided through the GitHub repository: [URaster Example Data on GitHub](https://github.com/changliao1025/uraster_data)


## 📚 Documentation

- [API Reference](https://uraster.readthedocs.io/en/latest/api.html)

## 📊 Supported Formats

- **Mesh formats**: GeoJSON, Shapefile, any OGR-supported vector format
- **Raster formats**: GeoTIFF, NetCDF, HDF5, any GDAL-supported raster format
- **Output formats**: Vectors (with computed statistics), PNG/JPG (visualizations), MP4/GIF (animations)

## 🙏 Acknowledgments

The model described in this repository was supported by the following:

* the U.S. Department of Energy Office of Science Biological and Environmental Research through the Earth System Development program as part of the Energy Exascale Earth System Model (E3SM) project.

* the Earth System Model Development and Regional and Global Model Analysis program areas of the U.S. Department of Energy, Office of Science, Biological and Environmental Research program as part of the multi-program, collaborative Integrated Coastal Modeling (ICoM) project.

* the Earth System Model Development and Regional and Global Model Analysis program areas of the U.S. Department of Energy, Office of Science, Biological and Environmental Research program as part of the multi-program, collaborative Interdisciplinary Research for Arctic Coastal Environments (InteRFACE) project.

A portion of this research was performed using PNNL Research Computing at Pacific Northwest National Laboratory.

PNNL is operated for DOE by Battelle Memorial Institute under contract DE-AC05-76RL01830.

## 🤝 Contributing & License

We welcome contributions! Please open an issue or submit a pull request on the GitHub repository.

**uraster** is distributed under the BSD 3-Clause License.

