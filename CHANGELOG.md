# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Initial package structure for PyPI release
- Comprehensive documentation and examples
- GitHub Actions CI/CD workflow
- Pre-commit hooks configuration

### Changed
- Updated README.md with detailed installation and usage instructions
- Fixed license consistency (BSD-3-Clause)
- Improved dependency specifications

### Fixed
- Added missing dependencies (pyearth, imageio)
- Corrected package metadata

## [0.1.0] - 2025-10-01

### Added
- Initial release of uraster package
- Core functionality for unstructured raster processing
- Support for zonal statistics on mesh geometries
- GDAL-native vector handling
- Interactive GeoVista 3D visualization
- Support for multiple raster and mesh formats
- Projection-aware operations
- Animation creation capabilities

### Features
- Structured raster to unstructured mesh conversion
- Zonal statistics computation
- 3D globe visualization with GeoVista
- Support for International Date Line crossing polygons
- Comprehensive error handling and crash detection
- Multiple remapping methods (nearest neighbor, weighted average)

### Dependencies
- numpy>=1.19.0
- gdal>=3.0.0
- pyearth>=0.1.0
- psutil>=5.0.0
- geovista>=0.3.0
- vtk==9.3.0
- imageio[ffmpeg]>=2.0.0