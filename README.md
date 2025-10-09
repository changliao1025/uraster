# uraster: Hybrid Geospatial Zonal Analysis and Visualization

Overview
uraster is the essential Python package for advanced geospatial data processing, designed to bridge the gap between structured raster data and vector-defined zones. It leverages GDAL/OGR for robust data handling.

The package provides a seamless, two-pronged workflow:

Analytical Mode (Zonal Statistics): Accurately calculates statistics (mean, min, max) from a raster (GeoTIFF, NetCDF) using GDAL and outputs the results as a standard Python object (e.g., list of dictionaries) for each polygon in the input vector.

Visualization Mode (Interactive 3D): Provides a direct API to generate a PyVista Unstructured Grid object, allowing for immediate, interactive, and georeferenced plotting of the statistical fields on a sphere using geovista.

✨ Core Features
GDAL-Native Vector Handling: Uses the robust GDAL/OGR engine for defining zones, reading vector files, and performing projection-aware operations.

Dual Output Workflow: Generates both a statistics-rich data structure (for analytical integration) and a visualization mesh (for 3D display).

Vector-Defined Zones: The input vector polygon defines the exact boundaries for both the statistical zones and the visualization mesh topology.

Projection-Aware Zonal Stats: Handles map projection differences to ensure accurate aggregation of raster values within each polygon.

Interactive GeoVista API: Offers simple functions to visualize the input raster and the output zonal statistics layer on a 3D sphere.

💻 Installation
uraster requires GDAL for vector handling and GeoVista (which relies on PyVista/VTK) for 3D visualization.

⚠️ GDAL Note: Installing GDAL's Python bindings can be complex via pip due to platform dependencies. We strongly recommend using Conda for a stable installation of GDAL.

# Install the core analysis package
pip install uraster

# Recommended: Install GDAL via Conda
conda install gdal

# Install visualization dependency
pip install geovista


🚀 Quick Start Example 1: Analytical Output (Zonal Statistics)
This example executes the zonal statistics and outputs a standard Python list containing the calculated statistics for each input polygon.




⚙️ Core API Methods (Summarized)
Method

Output Type

Description

uraster.zonal_stats_execute(...)

list[dict]

Analytical Output: Performs zonal statistics using GDAL and returns a list of dictionaries, where each dictionary contains the calculated statistics for one polygon.

uraster.to_geovista_mesh(...)

pyvista.UnstructuredGrid

Visualization Output: Converts the vector geometry and its attributes (from the analytical step) into a 3D mesh for interactive plotting on a sphere with geovista.

🤝 Contributing & License
We welcome contributions! Please open an issue or submit a pull request on the GitHub repository.

uraster is distributed under the MIT License.

