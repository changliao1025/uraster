# uraster: Structured Raster to Unstructured Mesh

## Overview

**uraster** is a Python package to convert or transfer structured raster dataset into unstructured mesh formats, designed to bridge the gap between structured raster data and unstructured mesh-based numerical models. It leverages GDAL/OGR for robust data handling.


[quickstart documentation](https://uraster.readthedocs.io/en/latest/quickstart.html)

## ✨ Core Features

- **GDAL-Native Vector Handling**: Uses the standard GDAL/OGR engine for defining unstructured mesh cells, and performing projection-aware geospatial operations.

- **Standard Vector I/O**: Instead of directly operating on various mesh standards, it utilizes standard geographic information system vector formats (e.g., GeoJSON) for mesh operations, ensuring broad compatibility. It supports transformation APIs between existing meshes and standard vector formats.

- **Projection-Aware Operations**: Handles (raster dateaset) map projection differences to ensure accurate aggregation of raster values within each polygon.

- **Interactive GeoVista API**: Offers simple functions to visualize the input and the output vector layers on a 3D sphere.

## 💻 Installation

uraster requires GDAL for vector handling and GeoVista (which relies on PyVista/VTK) for 3D visualization.

> ⚠️ **GDAL Note**: Installing GDAL's Python bindings can be complex via pip due to platform dependencies. We strongly recommend using Conda for a stable installation of GDAL and all dependencies.

### Install via Conda (Recommended)
```bash
# Create a new conda environment (recommended)
conda create -n uraster-env python=3.9
conda activate uraster-env

# Install uraster and all dependencies via conda
conda install -c conda-forge uraster
```

### Development Installation
```bash
git clone https://github.com/changliao1025/uraster.git
cd uraster
conda install -c conda-forge gdal geovista vtk=9.3.0 pyearth
conda develop .
```

## 🚀 Quick Start

### Example 1: Basic Zonal Statistics
This example demonstrates how to perform zonal statistics on unstructured mesh data:

```python
import uraster
from uraster.classes.uraster import uraster

# Configuration
config = {
    'sFilename_source_mesh': 'path/to/your/mesh.geojson',
    'aFilename_source_raster': ['path/to/your/raster.tif'],
    'sFilename_target_mesh': 'path/to/output/mesh_with_stats.geojson'
}

# Create uraster instance
processor = uraster(config)

# Setup and validate inputs
processor.setup(iFlag_verbose=True)

# Print input information
processor.report_inputs()

# Run zonal statistics
processor.run_remap(iFlag_verbose=True)

# Visualize results
processor.visualize_target_mesh(
    sVariable_in='mean',
    sFilename_out='output_visualization.png',
    sColormap='viridis'
)
```

### Example 2: Global Analysis with Animation
```python
import uraster
from uraster.classes.uraster import uraster

# Configuration for global analysis
config = {
    'sFilename_source_mesh': 'global_mesh.geojson',
    'aFilename_source_raster': ['global_dem.tif'],
    'sFilename_target_mesh': 'global_mesh_with_elevation.geojson'
}

processor = uraster(config)
processor.setup(iFlag_verbose=True)
processor.run_remap()

# Create rotating animation
processor.visualize_target_mesh(
    sVariable_in='mean',
    sFilename_out='global_elevation.mp4',
    sColormap='terrain',
    iFlag_create_animation=True,
    iAnimation_frames=360,
    sAnimation_format='mp4'
)
```

## 📊 Supported Formats

- **Mesh formats**: GeoJSON, Shapefile, any OGR-supported vector format
- **Raster formats**: GeoTIFF, NetCDF, HDF5, any GDAL-supported raster format
- **Output formats**: GeoJSON (with computed statistics), PNG/JPG (visualizations), MP4/GIF (animations)

## 🤝 Contributing & License

We welcome contributions! Please open an issue or submit a pull request on the GitHub repository.

**uraster** is distributed under the BSD 3-Clause License.

