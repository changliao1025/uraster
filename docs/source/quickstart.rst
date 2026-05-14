#####################
Quickstart
#####################

Installing and running uraster requires some basic knowledge of the Python ecosystem and Geographic Information System (GIS) concepts.

Follow these steps to quickly get started with uraster:

**1. Installation**

Install uraster using conda (recommended):

.. code-block:: bash

   # Create a new conda environment with uraster and all dependencies
   conda create -n uraster_env uraster vtk=9.3.0 -c conda-forge
   conda activate uraster_env

Or install from source for development:

.. code-block:: bash

   # Create and activate a new conda environment with dependencies
   conda create -n uraster_dev python=3.10 gdal numpy pyearth vtk=9.3.0 -c conda-forge
   conda activate uraster_dev

   # Clone and install uraster in development mode
   git clone https://github.com/changliao1025/uraster.git
   cd uraster
   pip install -e . --no-deps

For more details, see the :doc:`installation` guide.

**2. Basic Usage**

Here's a simple example to convert raster data to an unstructured mesh:

.. code-block:: python

   import uraster
   from uraster.classes.uraster import uraster

   # Configuration dictionary
   config = {
       'sFilename_source_mesh': 'path/to/your/mesh.geojson',
       'aFilename_source_raster': ['path/to/your/raster.tif'],
       'sFilename_target_mesh': 'path/to/output/mesh_with_stats.geojson',
       'iFlag_remap_method': 1,  # 1=nearest neighbor, 2=bilinear, 3=weighted average
       'sField_unique_id': 'cellid',  # Field name for unique cell IDs
       'iFlag_global': 1,  # 1 for global datasets, 0 for regional
       'iFlag_polar': 0   # 1 if data includes polar regions, 0 otherwise
   }

   # Create uraster instance
   processor = uraster(config)

   # Setup and validate inputs
   processor.setup(iFlag_verbose_in=True)

   # Run zonal statistics to transfer raster values to mesh
   processor.run_remap(iFlag_verbose_in=True)

**2a. Advanced Configuration**

For more control over processing:

.. code-block:: python

   # Advanced configuration with discrete data support
   config_advanced = {
       'sFilename_source_mesh': 'path/to/your/mesh.geojson',
       'aFilename_source_raster': ['path/to/discrete_raster.tif'],
       'sFilename_target_mesh': 'path/to/output/mesh_with_stats.geojson',
       'iFlag_remap_method': 1,  # Only method 1 supported for discrete data
       'sField_unique_id': 'cellid'
   }

   processor = uraster(config_advanced)
   processor.setup(iFlag_verbose_in=True)

   # Run with discrete data handling
   processor.run_remap(
       iFlag_discrete_in=True,  # Enable discrete data mode
       iFlag_stat_in=False,     # Disable statistics for discrete data
       iFlag_verbose_in=True
   )

**2b. Memory Management**

uraster includes intelligent memory management for large datasets:

.. code-block:: python

   # Configure memory management
   processor.set_cache_threshold(threshold_mb=200)  # Set cache threshold to 200MB

   # Get cache information
   cache_info = processor._get_cache_info()
   print(f"Cached files: {cache_info['cached_files']}")

   # Use context manager for data-intensive operations
   with processor.get_sraster_with_data('large_raster.tif') as raster:
       raster.read_data()  # Process large raster data
       # Automatic cleanup when exiting context

**3. Visualization (Optional)**

Visualize your results using the enhanced GeoVista 3D visualization:

.. code-block:: python

   # Basic visualization of the output mesh with statistics
   processor.visualize_target_mesh(
       sVariable_in='mean',
       sFilename_out='output_visualization.png',
       sColormap='viridis'
   )

   # Advanced visualization with custom focus and styling
   processor.visualize_target_mesh(
       sVariable_in='mean',
       sUnit_in='mm/year',  # Unit label for colorbar
       sFilename_out='advanced_viz.png',
       dLongitude_focus_in=-100.0,  # Focus on Americas
       dLatitude_focus_in=45.0,
       dZoom_factor=0.8,
       iFlag_show_coastlines=True,
       iFlag_show_graticule=True,
       sColormap='plasma'
   )

   # Create rotating animation
   processor.visualize_target_mesh(
       sVariable_in='mean',
       sFilename_out='rotation_animation.mp4',
       iFlag_create_animation=True,
       iAnimation_frames=72,  # 5° per frame for smooth rotation
       sAnimation_format='mp4',
       iFlag_verbose_in=True
   )

   # Visualize source mesh topology
   processor.visualize_source_mesh(
       sFilename_out='source_mesh.png',
       dLongitude_focus_in=0.0,
       dLatitude_focus_in=0.0,
       iFlag_show_coastlines=True
   )

**4. Examples**

Check the `examples/` directory for more detailed usage scenarios:

- Global examples: Working with global datasets and MPAS/DGGS meshes
- Regional examples: Working with TIN (Triangulated Irregular Network) meshes

If you encounter any issues, refer to the FAQ or submit a GitHub issue (https://github.com/changliao1025/uraster/issues).