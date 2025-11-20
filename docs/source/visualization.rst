##############
Visualization
##############

uraster provides powerful 3D visualization capabilities using GeoVista for interactive globe rendering. This allows you to visualize both input mesh topology and output results with computed statistics.

Prerequisites
=============

To use visualization features, you need to install GeoVista:

.. code-block:: bash

   pip install geovista

For animation creation, also install imageio:

.. code-block:: bash

   pip install imageio[ffmpeg]

Basic Visualization
==================

Source Mesh Visualization
-------------------------

Visualize the topology of your input mesh to verify geometry and coverage:

.. code-block:: python

   from uraster.classes.uraster import uraster

   config = {
       'sFilename_source_mesh': 'path/to/your/mesh.geojson',
       'aFilename_source_raster': ['path/to/your/raster.tif']
   }

   processor = uraster(config)
   processor.setup(iFlag_verbose_in=True)

   # Interactive visualization (opens viewer window)
   processor.visualize_source_mesh()

   # Save static image
   processor.visualize_source_mesh(
       sFilename_out='source_mesh.png',
       dLongitude_focus_in=0.0,  # Focus longitude
       dLatitude_focus_in=0.0,   # Focus latitude
       dZoom_factor=0.7,         # Zoom level
       iFlag_show_coastlines=True,
       iFlag_show_graticule=True
   )

Target Mesh Visualization with Data
----------------------------------

Visualize computed statistics on the mesh after processing:

.. code-block:: python

   # First run the remapping to generate statistics
   processor.run_remap(iFlag_verbose_in=True)

   # Visualize the results
   processor.visualize_target_mesh(
       sVariable_in='mean',  # Variable to visualize
       sUnit_in='mm',        # Unit for colorbar
       sFilename_out='results.png',
       sColormap='viridis'   # Matplotlib colormap
   )

Advanced Visualization Options
=============================

Custom Focus and Styling
------------------------

Control camera position, zoom, and visual elements:

.. code-block:: python

   processor.visualize_target_mesh(
       sVariable_in='mean',
       sUnit_in='kg/m²',
       sFilename_out='focused_view.png',
       dLongitude_focus_in=-100.0,  # Focus on Americas
       dLatitude_focus_in=45.0,     # Focus on mid-latitudes
       dZoom_factor=0.8,            # Closer zoom
       iFlag_show_coastlines=True,  # Show coastlines
       iFlag_show_graticule=True,   # Show coordinate grid
       sColormap='plasma',          # Different colormap
       iFlag_verbose_in=True
   )

Available Colormaps
------------------

uraster supports all matplotlib colormaps. Popular choices include:

- ``'viridis'`` - Default, perceptually uniform
- ``'plasma'`` - Purple-pink-yellow
- ``'coolwarm'`` - Blue-white-red diverging
- ``'RdYlBu'`` - Red-yellow-blue diverging
- ``'jet'`` - Rainbow (not recommended for scientific data)
- ``'terrain'`` - Earth-like colors for topography

Animation Creation
=================

Creating Rotating Animations
----------------------------

Generate smooth rotating animations to showcase global datasets:

.. code-block:: python

   # Create a 360° rotation animation
   processor.visualize_target_mesh(
       sVariable_in='mean',
       sUnit_in='mm/year',
       sFilename_out='global_rotation.mp4',
       dLongitude_focus_in=0.0,     # Starting longitude
       dLatitude_focus_in=0.0,      # Camera latitude
       dZoom_factor=0.75,           # Globe zoom level
       iFlag_create_animation=True, # Enable animation
       iAnimation_frames=72,        # 72 frames = 5° per frame
       dAnimation_speed=1.0,        # Animation speed multiplier
       sAnimation_format='mp4',     # Output format
       iFlag_show_coastlines=True,
       iFlag_show_graticule=True,
       sColormap='viridis',
       iFlag_verbose_in=True
   )

Animation Parameters
-------------------

.. list-table:: Animation Parameters
   :widths: 25 15 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``iFlag_create_animation``
     - bool
     - Enable animation mode
   * - ``iAnimation_frames``
     - int
     - Number of frames for 360° rotation (default: 360)
   * - ``dAnimation_speed``
     - float
     - Speed multiplier for rotation (default: 1.0)
   * - ``sAnimation_format``
     - str
     - Output format: 'mp4', 'gif', 'avi' (default: 'mp4')

Visualization Parameters Reference
=================================

Common Parameters
----------------

.. list-table:: Common Visualization Parameters
   :widths: 25 15 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``sFilename_out``
     - str
     - Output file path (None for interactive)
   * - ``dLongitude_focus_in``
     - float
     - Camera focus longitude (-180 to 180)
   * - ``dLatitude_focus_in``
     - float
     - Camera focus latitude (-90 to 90)
   * - ``dZoom_factor``
     - float
     - Zoom level (higher = closer)
   * - ``iFlag_show_coastlines``
     - bool
     - Show coastline overlay
   * - ``iFlag_show_graticule``
     - bool
     - Show coordinate grid
   * - ``iFlag_verbose_in``
     - bool
     - Enable verbose logging

Target Mesh Specific Parameters
------------------------------

.. list-table:: Target Mesh Visualization Parameters
   :widths: 25 15 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``sVariable_in``
     - str
     - Variable field to visualize ('mean', 'min', 'max', 'std')
   * - ``sUnit_in``
     - str
     - Unit label for colorbar
   * - ``sColormap``
     - str
     - Matplotlib colormap name

Best Practices
=============

Performance Tips
---------------

1. **Large Datasets**: Use static image output instead of interactive mode for very large meshes
2. **Animation Quality**: Use 36-72 frames for smooth animations (10°-5° per frame)
3. **Memory Management**: Close interactive viewers when done to free GPU memory

Visual Quality Tips
------------------

1. **Colormap Selection**: Use perceptually uniform colormaps like 'viridis' for scientific data
2. **Focus Selection**: Choose focus points that highlight important features in your data
3. **Zoom Levels**: Use zoom factors between 0.5-1.0 for global views, higher for regional focus
4. **Coastlines**: Enable coastlines for geographic context, especially for oceanographic data

Output Formats
=============

Supported Formats
----------------

**Static Images:**
- PNG (recommended for publications)
- JPG (smaller file size)
- SVG (vector format, scalable)

**Animations:**
- MP4 (recommended, widely compatible)
- GIF (larger file size, universal support)
- AVI (high quality, larger files)

File Size Considerations
-----------------------

- PNG: Best quality, moderate size
- JPG: Smaller size, some quality loss
- MP4: Efficient compression for animations
- GIF: Large files but universal browser support

Troubleshooting
==============

Common Issues
------------

**"GeoVista not available"**
   Install GeoVista: ``pip install geovista``

**"No display available"**
   Use ``sFilename_out`` parameter to save files instead of interactive display

**Animation creation fails**
   Install imageio with ffmpeg: ``pip install imageio[ffmpeg]``

**GPU/Memory issues**
   - Reduce mesh complexity for visualization
   - Use static output instead of interactive mode
   - Close visualization windows promptly

**Empty visualization**
   - Check that target mesh file exists and contains data
   - Verify variable name exists in the mesh file
   - Ensure mesh topology was built successfully

Performance Optimization
-----------------------

For large datasets:

1. Use appropriate zoom levels to focus on areas of interest
2. Consider mesh simplification for visualization purposes
3. Use static output formats for batch processing
4. Monitor memory usage with large global meshes

Examples
========

Regional Focus Example
---------------------

.. code-block:: python

   # Focus on North America
   processor.visualize_target_mesh(
       sVariable_in='precipitation',
       sUnit_in='mm/day',
       sFilename_out='north_america_precip.png',
       dLongitude_focus_in=-100.0,
       dLatitude_focus_in=45.0,
       dZoom_factor=1.2,
       sColormap='Blues',
       iFlag_show_coastlines=True
   )

Global Animation Example
-----------------------

.. code-block:: python

   # Create smooth global rotation
   processor.visualize_target_mesh(
       sVariable_in='temperature',
       sUnit_in='°C',
       sFilename_out='global_temp_rotation.mp4',
       iFlag_create_animation=True,
       iAnimation_frames=60,  # 6° per frame
       sAnimation_format='mp4',
       sColormap='RdYlBu_r',  # Reversed red-yellow-blue
       iFlag_show_coastlines=True,
       iFlag_verbose_in=True
   )