uraster: Structured Raster to Unstructured Mesh
=================================================

.. image:: https://img.shields.io/pypi/v/uraster.svg
   :target: https://pypi.org/project/uraster/
   :alt: PyPI version

.. image:: https://img.shields.io/conda/vn/conda-forge/uraster.svg
   :target: https://anaconda.org/conda-forge/uraster
   :alt: Conda version

.. image:: https://github.com/changliao1025/uraster/workflows/CI/badge.svg
   :target: https://github.com/changliao1025/uraster/actions
   :alt: CI Status

.. image:: https://codecov.io/gh/changliao1025/uraster/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/changliao1025/uraster
   :alt: Coverage Status

**uraster** is a Python package to convert or transfer structured raster dataset into unstructured mesh formats, designed to bridge the gap between structured raster data and unstructured mesh-based numerical models. It leverages GDAL/OGR for robust data handling.

Core Features
-------------

- **GDAL-Native Vector Handling**: Uses the standard GDAL/OGR engine for defining unstructured mesh cells, and performing projection-aware geospatial operations.

- **Standard Vector I/O**: Instead of directly operating on various mesh standards, it utilizes standard geographic information system vector formats (e.g., GeoJSON) for mesh operations, ensuring broad compatibility. It supports transformation APIs between existing meshes and standard vector formats.

- **Projection-Aware Operations**: Handles (raster dataset) map projection differences to ensure accurate aggregation of raster values within each polygon.

- **Interactive GeoVista API**: Offers simple functions to visualize the input and the output vector layers on a 3D sphere.

Quick Start
-----------

.. code-block:: python

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

   # Run zonal statistics
   processor.run_remap(iFlag_verbose=True)

   # Visualize results
   processor.visualize_target_mesh(
       sVariable_in='mean',
       sFilename_out='output_visualization.png',
       sColormap='viridis'
   )

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   user_guide
   api_reference
   examples
   contributing
   changelog

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`