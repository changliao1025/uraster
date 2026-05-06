# Installation

## From PyPI (when published)

```bash
pip install uraster
```

## From Source

```bash
git clone https://github.com/changliao1025/uraster.git
cd uraster
pip install -e .
```

## Dependencies

### Core Dependencies
- Python >= 3.8
- numpy >= 1.19.0
- GDAL >= 3.0.0
- osgeo (GDAL Python bindings)
- pyearth >= 0.2.1

### Optional Dependencies
- pyearth >= 0.1.1 (for Earth data handling)
- geovista (for 3D visualization)
- imageio[ffmpeg] (for animation creation)

### Installing Optional Dependencies

For full functionality including 3D visualization and animations:

```bash
pip install geovista imageio[ffmpeg]
```

Or install with conda:

```bash
conda install -c conda-forge geovista imageio
```

### Installing GDAL

GDAL can be tricky to install. Here are platform-specific instructions:

**Ubuntu/Debian:**
```bash
sudo apt-get install gdal-bin libgdal-dev
pip install gdal
```

**macOS (with Homebrew):**
```bash
brew install gdal
pip install gdal
```

**Windows:**
```bash
# Use conda for easier GDAL installation
conda install -c conda-forge gdal
```