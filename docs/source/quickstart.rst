#####################
Quickstart
#####################

Installing and running uraster requires some basic knowledge of the Python ecosystem and Geographic Information System (GIS) concepts.

Follow these steps to quickly get started with uraster:

**1. Installation**

Install uraster using pip:

.. code-block:: bash

   pip install uraster

Or install from source:

.. code-block:: bash

   git clone https://github.com/changliao1025/uraster.git
   cd uraster
   pip install -e .

**2. Basic Usage**

Here's a simple example to convert raster data to an unstructured mesh:

.. code-block:: python

   import uraster
   from uraster.classes.uraster import uraster

   # Configuration dictionary
   config = {
       'sFilename_source_mesh': 'path/to/your/mesh.geojson',
       'aFilename_source_raster': ['path/to/your/raster.tif'],
       'sFilename_target_mesh': 'path/to/output/mesh_with_stats.geojson'
   }

   # Create uraster instance
   processor = uraster(config)

   # Setup and validate inputs
   processor.setup(iFlag_verbose=True)

   # Run zonal statistics to transfer raster values to mesh
   processor.run_remap(iFlag_verbose=True)

**3. Visualization (Optional)**

Visualize your results:

.. code-block:: python

   # Visualize the output mesh with statistics
   processor.visualize_target_mesh(
       sVariable_in='mean',
       sFilename_out='output_visualization.png',
       sColormap='viridis'
   )

**4. Examples**

Check the `examples/` directory for more detailed usage scenarios:

- Global examples: Working with global datasets and MPAS/DGGS meshes
- Regional examples: Working with TIN (Triangulated Irregular Network) meshes

If you encounter any issues, refer to the FAQ or submit a GitHub issue (https://github.com/changliao1025/uraster/issues).