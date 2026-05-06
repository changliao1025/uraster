# Installation

URaster requires GDAL for vector handling, which can be complex to install via pip due to platform dependencies. We strongly recommend using Conda for a stable installation of GDAL and all dependencies.

## Install via Conda (Recommended)

```bash
# Create a new conda environment with uraster and all dependencies
conda create -n uraster_env uraster vtk=9.3.0 -c conda-forge
conda activate uraster_env
```

This single command will automatically install uraster and all its dependencies, including GDAL, numpy, pyearth, and VTK 9.3.0 for visualization compatibility.

## Install from Source with Conda

If you want to install from source for development:

```bash
# Create and activate a new conda environment with dependencies
conda create -n uraster_dev python=3.10 gdal numpy pyearth vtk=9.3.0 -c conda-forge
conda activate uraster_dev

# Clone and install uraster in development mode
git clone https://github.com/changliao1025/uraster.git
cd uraster
pip install -e . --no-deps
```

## Dependencies

### Core Dependencies
- Python >= 3.8
- numpy >= 1.19.0
- GDAL >= 3.0.0
- osgeo (GDAL Python bindings)
- pyearth >= 0.2.1

### Optional Dependencies
- geovista (for 3D visualization)
- imageio (for animation creation)
- vtk = 9.3.0 (for PyVista compatibility)

### Installing Optional Dependencies

For full functionality including 3D visualization and animations:

```bash
conda install -c conda-forge geovista imageio
```

> **Note**: Due to compatibility issues in PyVista, we recommend using VTK version 9.3.0 for stable visualization functionality.
