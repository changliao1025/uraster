# uraster: Structured Raster to Unstructured Mesh

## Overview

**uraster** is a Python package to convert or transfer structured raster dataset into unstructured mesh formats, designed to bridge the gap between structured raster data and unstructured mesh-based numerical models. It leverages GDAL/OGR for robust data handling.


## ✨ Core Features

- **GDAL-Native Vector Handling**: Uses the standard GDAL/OGR engine for defining unstructured mesh cells, and performing projection-aware geospatial operations.

- **Standard Vector I/O**: Instead of directly operating on various mesh standards, it utilizes standard geographic information system vector formats (e.g., GeoJSON) for mesh operations, ensuring broad compatibility. It supports transformation APIs between existing meshes and standard vector formats.

- **Projection-Aware Operations**: Handles (raster dateaset) map projection differences to ensure accurate aggregation of raster values within each polygon.

- **Interactive GeoVista API**: Offers simple functions to visualize the input and the output vector layers on a 3D sphere.

## 💻 Installation

uraster requires GDAL for vector handling and GeoVista (which relies on PyVista/VTK) for 3D visualization.

> ⚠️ **GDAL Note**: Installing GDAL's Python bindings can be complex via pip due to platform dependencies. We strongly recommend using Conda for a stable installation of GDAL.

### Install uraster via Conda (Recommended)
```bash
# Recommended: Install uraster via Conda
conda install uraster
```

## 🚀 Quick Start

### Example 1: Analytical Output (Zonal Statistics)
This example executes the zonal statistics and outputs a standard Python list containing the calculated statistics for each input polygon.

```python
import uraster

# TODO: Add example code here
```

## 🤝 Contributing & License

We welcome contributions! Please open an issue or submit a pull request on the GitHub repository.

**uraster** is distributed under the MIT License.

