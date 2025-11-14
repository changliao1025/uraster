"""
Visualization module for uraster class.

This module contains all visualization-related methods that were moved from the main uraster class
to reduce the size of the main uraster.py file and improve code organization.

Features:
- 3D mesh visualization using GeoVista
- Interactive and static rendering modes
- Animation support with rotation and camera movement
- Comprehensive error handling and validation
- Support for multiple output formats
"""

import os
import logging
import traceback
import math
from typing import Optional, List, Tuple, Union, Dict, Any
import numpy as np
from osgeo import gdal, ogr
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
from pyearth.gis.geometry.extract_unique_vertices_and_connectivity import extract_unique_vertices_and_connectivity
from uraster.classes.sraster import sraster

gdal.UseExceptions()

# Set up logging
logger = logging.getLogger(__name__)
CRS = "EPSG:4326"

# Constants for visualization
DEFAULT_EARTH_RADIUS = 1.0
DEFAULT_CAMERA_DISTANCE_MULTIPLIER = 3.0
DEFAULT_ZOOM_FACTOR = 0.7
VALID_ANIMATION_FORMATS = ['mp4', 'gif', 'avi']
VALID_IMAGE_FORMATS = ['.png', '.jpg', '.jpeg', '.svg', '.tif', '.tiff']
COORDINATE_BOUNDS = {'longitude': (-180, 180), 'latitude': (-90, 90)}

class VisualizationConfig:
    """Configuration class for visualization parameters."""

    def __init__(self,
                 longitude_focus: float = 0.0,
                 latitude_focus: float = 0.0,
                 zoom_factor: float = DEFAULT_ZOOM_FACTOR,
                 show_coastlines: bool = True,
                 show_graticule: bool = True,
                 colormap: str = 'viridis',
                 coastline_color: str = 'black',
                 coastline_width: float = 1.0,
                 verbose: bool = False):
        self.longitude_focus = self._validate_longitude(longitude_focus)
        self.latitude_focus = self._validate_latitude(latitude_focus)
        self.zoom_factor = self._validate_zoom_factor(zoom_factor)
        self.show_coastlines = show_coastlines
        self.show_graticule = show_graticule
        self.colormap = colormap
        self.coastline_color = coastline_color
        self.coastline_width = coastline_width
        self.verbose = verbose

    def _validate_longitude(self, lon: float) -> float:
        """Validate and clamp longitude to valid range."""
        if not (-180 <= lon <= 180):
            logger.warning(f'Longitude {lon} out of range [-180, 180], clamping')
            return np.clip(lon, -180, 180)
        return lon

    def _validate_latitude(self, lat: float) -> float:
        """Validate and clamp latitude to valid range."""
        if not (-90 <= lat <= 90):
            logger.warning(f'Latitude {lat} out of range [-90, 90], clamping')
            return np.clip(lat, -90, 90)
        return lat

    def _validate_zoom_factor(self, zoom: float) -> float:
        """Validate zoom factor."""
        if zoom <= 0:
            logger.warning(f'Invalid zoom factor {zoom}, using default {DEFAULT_ZOOM_FACTOR}')
            return DEFAULT_ZOOM_FACTOR
        return zoom

class CameraController:
    """Handles camera positioning and movement calculations."""

    @staticmethod
    def calculate_camera_position(longitude: float, latitude: float,
                                zoom_factor: float = DEFAULT_ZOOM_FACTOR) -> Tuple[List[float], List[float]]:
        """
        Calculate camera position and focal point from geographic coordinates.

        Args:
            longitude: Longitude in degrees
            latitude: Latitude in degrees
            zoom_factor: Camera zoom level

        Returns:
            Tuple of (focal_point, camera_position) as [x, y, z] lists
        """
        # Convert to radians
        lon_rad = math.radians(longitude)
        lat_rad = math.radians(latitude)

        # Calculate positions
        earth_radius = DEFAULT_EARTH_RADIUS
        camera_distance = earth_radius * DEFAULT_CAMERA_DISTANCE_MULTIPLIER

        # Focal point on Earth surface
        x_focal = earth_radius * math.cos(lat_rad) * math.cos(lon_rad)
        y_focal = earth_radius * math.cos(lat_rad) * math.sin(lon_rad)
        z_focal = earth_radius * math.sin(lat_rad)

        # Camera position away from Earth
        x_camera = camera_distance * math.cos(lat_rad) * math.cos(lon_rad)
        y_camera = camera_distance * math.cos(lat_rad) * math.sin(lon_rad)
        z_camera = camera_distance * math.sin(lat_rad)

        focal_point = [x_focal, y_focal, z_focal]
        camera_position = [x_camera, y_camera, z_camera]

        return focal_point, camera_position

    @staticmethod
    def validate_camera_setup(focal_point: List[float], camera_position: List[float]) -> bool:
        """Validate camera setup to ensure proper positioning."""
        distance = math.sqrt(
            sum((c - f)**2 for c, f in zip(camera_position, focal_point))
        )
        return distance >= 0.1  # Minimum distance threshold

class AnimationConfig:
    """Configuration class for animation parameters."""

    def __init__(self,
                 frames: int = 36,
                 speed: float = 1.0,
                 format: str = 'mp4',
                 amplitude_deg: float = 20.0,
                 cycles: float = 1.0,
                 phase: float = 0.0):
        # Convert to int if string and validate
        try:
            frames_int = int(frames) if isinstance(frames, str) else frames
            self.frames = max(1, frames_int)  # Ensure at least 1 frame
        except (ValueError, TypeError):
            logger.warning(f'Invalid frames value {frames}, using default 36')
            self.frames = 36

        # Convert to float if string and validate
        try:
            speed_float = float(speed) if isinstance(speed, str) else speed
            self.speed = max(0.1, speed_float)  # Ensure positive speed
        except (ValueError, TypeError):
            logger.warning(f'Invalid speed value {speed}, using default 1.0')
            self.speed = 1.0

        self.format = format.lower()
        self.amplitude_deg = amplitude_deg
        self.cycles = cycles
        self.phase = phase

        # Validate format
        if self.format not in VALID_ANIMATION_FORMATS:
            logger.warning(f'Invalid animation format {format}, using mp4')
            self.format = 'mp4'

def _validate_mesh_data(uraster_instance) -> bool:
    """
    Validate that mesh data is available and properly formatted.

    Args:
        uraster_instance: The uraster instance to validate

    Returns:
        bool: True if mesh data is valid, False otherwise
    """
    if uraster_instance.aVertex_longititude is None or uraster_instance.aVertex_latitude is None:
        logger.error('Mesh vertices not available. Build mesh topology first.')
        return False

    if uraster_instance.aConnectivity is None:
        logger.error('Mesh connectivity not available. Build mesh topology first.')
        return False

    if len(uraster_instance.aVertex_longititude) == 0 or len(uraster_instance.aVertex_latitude) == 0:
        logger.error('Mesh vertices are empty.')
        return False

    if uraster_instance.aConnectivity.size == 0:
        logger.error('Mesh connectivity is empty.')
        return False

    return True

def _validate_output_path(sFilename: Optional[str]) -> bool:
    """
    Validate output file path and create directories if needed.

    Args:
        sFilename: Output file path to validate

    Returns:
        bool: True if path is valid and accessible, False otherwise
    """
    if sFilename is None:
        return True  # Interactive mode, no file validation needed

    if not isinstance(sFilename, str) or not sFilename.strip():
        logger.error('Output filename must be a non-empty string')
        return False

    # Check output directory exists
    sOutput_dir = os.path.dirname(sFilename)
    if sOutput_dir and not os.path.exists(sOutput_dir):
        try:
            os.makedirs(sOutput_dir, exist_ok=True)
            logger.info(f'Created output directory: {sOutput_dir}')
        except Exception as e:
            logger.error(f'Cannot create output directory {sOutput_dir}: {e}')
            return False

    # Check supported file extensions
    sFile_ext = os.path.splitext(sFilename)[1].lower()
    aAll_valid_extensions = VALID_IMAGE_FORMATS + [f'.{fmt}' for fmt in VALID_ANIMATION_FORMATS]
    if sFile_ext not in aAll_valid_extensions:
        logger.warning(f'File extension {sFile_ext} may not be supported. '
                      f'Recommended: {", ".join(VALID_IMAGE_FORMATS)}')

    return True

def _setup_geovista_plotter(iFlag_off_screen: bool = False, iFlag_verbose: bool = False):
    """
    Set up GeoVista plotter with error handling.

    Args:
        iFlag_off_screen: Whether to create off-screen plotter
        iFlag_verbose: Enable verbose logging

    Returns:
        GeoVista plotter instance or None if failed
    """
    try:
        import geovista as gv
        if iFlag_verbose:
            logger.info('GeoVista library imported successfully')
    except ImportError as e:
        logger.error('GeoVista library not available. Install with: pip install geovista')
        logger.error(f'Import error: {e}')
        return None

    try:
        if iFlag_off_screen:
            pPlotter = gv.GeoPlotter(off_screen=True)
            if iFlag_verbose:
                logger.debug('Created off-screen plotter')
        else:
            pPlotter = gv.GeoPlotter()
            if iFlag_verbose:
                logger.debug('Created interactive plotter')
        return pPlotter
    except Exception as e:
        logger.error(f'Failed to create GeoVista plotter: {e}')
        logger.error('This may be due to missing graphics context or display')
        if iFlag_off_screen:
            logger.error('For headless systems, ensure proper OpenGL/Mesa setup')
        else:
            logger.error('For interactive mode, ensure display environment (X11/Wayland) is available')
        return None

def _add_geographic_context(pPlotter, pConfig: VisualizationConfig):
    """
    Add geographic context (coastlines, graticule, axes) to plotter.

    Args:
        pPlotter: GeoVista plotter instance
        pConfig: Visualization configuration
    """
    # Add coastlines
    if pConfig.show_coastlines:
        try:
            # You can set coastline color using the 'color' parameter
            # Common options: 'black', 'white', 'red', 'blue', 'gray', etc.
            # You can also use RGB tuples like (1.0, 0.0, 0.0) for red
            pPlotter.add_coastlines(color=pConfig.coastline_color, line_width=pConfig.coastline_width)
            if pConfig.verbose:
                logger.debug(f'Added coastlines overlay (color: {pConfig.coastline_color}, width: {pConfig.coastline_width})')
        except Exception as e:
            logger.warning(f'Could not add coastlines: {e}')

    # Add coordinate axes
    try:
        pPlotter.add_axes()
        if pConfig.verbose:
            logger.debug('Added coordinate axes')
    except Exception as e:
        logger.warning(f'Could not add axes: {e}')

    # Add graticule (coordinate grid)
    if pConfig.show_graticule:
        try:
            pPlotter.add_graticule(show_labels=True)
            if pConfig.verbose:
                logger.debug('Added coordinate graticule with labels')
        except Exception as e:
            logger.warning(f'Could not add graticule: {e}')

def _configure_camera(pPlotter, pConfig: VisualizationConfig) -> bool:
    """
    Configure camera position and orientation.

    Args:
        pPlotter: GeoVista plotter instance
        pConfig: Visualization configuration

    Returns:
        bool: True if camera configured successfully, False otherwise
    """
    try:
        aFocal_point, aCamera_position = CameraController.calculate_camera_position(
            pConfig.longitude_focus, pConfig.latitude_focus, pConfig.zoom_factor
        )

        # Validate camera setup
        if not CameraController.validate_camera_setup(aFocal_point, aCamera_position):
            logger.warning('Camera and focal point are too close, using default view')
            raise ValueError('Invalid camera positioning')

        pPlotter.camera.focal_point = aFocal_point
        pPlotter.camera.position = aCamera_position
        pPlotter.camera.zoom(pConfig.zoom_factor)
        pPlotter.camera.up = [0, 0, 1]  # Z-up orientation

        if pConfig.verbose:
            logger.debug(f'Camera configured successfully:')
            logger.debug(f'  Focal point: {aFocal_point}')
            logger.debug(f'  Camera position: {aCamera_position}')

        return True

    except Exception as e:
        logger.warning(f'Error setting camera position: {e}. Using default view.')
        try:
            pPlotter.reset_camera()
        except Exception:
            pass
        return False

def visualize_source_mesh(self,
                          sFilename_out: Optional[str] = None,
                          dLongitude_focus_in: Optional[float] = 0.0,
                          dLatitude_focus_in: Optional[float] = 0.0,
                          dZoom_factor: float = 0.7,
                          iFlag_show_coastlines: bool = True,
                          iFlag_show_graticule: bool = True,
                          sCoastline_color: str = 'black',
                          dCoastline_width: float = 1.0,
                          iFlag_verbose: bool = False) -> bool:
    """
    Visualize the source mesh topology using GeoVista 3D globe rendering.

    Creates an interactive or saved 3D visualization of the unstructured mesh
    with proper geographic context including coastlines and coordinate grid.

    Args:
        sFilename_out: Output screenshot file path. If None, displays interactive viewer.
            Supports formats: .png, .jpg, .svg
        dLongitude_focus_in: Camera focal point longitude in degrees (-180 to 180).
            Default is 0.0 (prime meridian).
        dLatitude_focus_in: Camera focal point latitude in degrees (-90 to 90).
            Default is 0.0 (equator).
        dZoom_factor: Camera zoom level. Higher values zoom in. Default is 0.7.
        iFlag_show_coastlines: Show coastline overlay. Default is True.
        iFlag_show_graticule: Show coordinate grid with labels. Default is True.
        sCoastline_color: Color for coastlines. Default is 'black'.
            Examples: 'white', 'red', 'blue', 'gray', or RGB tuples like (1.0, 0.0, 0.0).
        dCoastline_width: Line width for coastlines. Default is 1.0.
        iFlag_verbose: If True, print detailed progress messages. Default is False.

    Returns:
        True if visualization successful, False otherwise

    Note:
        - Requires 'geovista' package: pip install geovista
        - Interactive mode requires display environment
        - Mesh topology must be built before visualization (call rebuild_mesh_topology first)
    """
    # Validate inputs using new utility functions
    if not _validate_mesh_data(self):
        return False

    if not _validate_output_path(sFilename_out):
        return False

    # Create configuration object
    config = VisualizationConfig(
        longitude_focus=dLongitude_focus_in,
        latitude_focus=dLatitude_focus_in,
        zoom_factor=dZoom_factor,
        show_coastlines=iFlag_show_coastlines,
        show_graticule=iFlag_show_graticule,
        coastline_color=sCoastline_color,
        coastline_width=dCoastline_width,
        verbose=iFlag_verbose
    )

    try:
        # Import and setup GeoVista
        import geovista as gv
        if config.verbose:
            logger.info('Creating mesh visualization...')
            logger.info(f'  - Vertices: {len(self.aVertex_longititude)}')
            logger.info(f'  - Connectivity shape: {self.aConnectivity.shape}')
            logger.info(f'  - Focus: ({config.longitude_focus:.2f}°, {config.latitude_focus:.2f}°)')
            logger.info(f'  - Zoom factor: {config.zoom_factor}')

        # Validate connectivity array structure
        if self.aConnectivity.ndim != 2:
            logger.error(f'Connectivity array must be 2D, got {self.aConnectivity.ndim}D')
            return False

        # Create masked connectivity array (mask invalid indices)
        connectivity_masked = np.ma.masked_where(
            self.aConnectivity == -1,
            self.aConnectivity
        )

        # Validate connectivity indices
        valid_connectivity = self.aConnectivity[self.aConnectivity >= 0]
        if len(valid_connectivity) > 0:
            max_vertex_idx = len(self.aVertex_longititude) - 1
            if np.max(valid_connectivity) > max_vertex_idx:
                logger.error(f'Connectivity contains invalid vertex index: '
                           f'max={np.max(valid_connectivity)}, vertices={len(self.aVertex_longititude)}')
                return False

        # Transform to GeoVista unstructured mesh
        mesh = gv.Transform.from_unstructured(
            self.aVertex_longititude,
            self.aVertex_latitude,
            connectivity=connectivity_masked,
            crs=CRS
        )

        # Validate cell data array length matches mesh cells
        if len(self.aCellID) != mesh.n_cells:
            logger.error(f'Cell ID array length ({len(self.aCellID)}) does not match '
                        f'mesh cells ({mesh.n_cells})')
            return False

        # Prepare mesh metadata
        name = 'Mesh Cell ID'
        mesh.cell_data[name] = self.aCellID

        if config.verbose:
            logger.info(f'Created GeoVista mesh with {mesh.n_cells} cells and {mesh.n_points} points')

        # Setup plotter
        pPlotter = _setup_geovista_plotter(iFlag_off_screen=(sFilename_out is not None), iFlag_verbose=config.verbose)
        if pPlotter is None:
            return False

        # Configure scalar bar (colorbar) appearance
        sargs = {
            "title": name,
            "shadow": True,
            "title_font_size": 10,
            "label_font_size": 10,
            "fmt": "%.0f",  # Integer formatting for cell IDs
            "n_labels": 5,
        }

        # Add mesh to plotter
        pPlotter.add_mesh(mesh, scalars=name, scalar_bar_args=sargs)

        # Configure camera
        _configure_camera(pPlotter, config)

        # Add geographic context
        _add_geographic_context(pPlotter, config)

        # Output or display
        return _handle_visualization_output(pPlotter, sFilename_out, config.verbose)

    except ImportError as e:
        logger.error('GeoVista library not available. Install with: pip install geovista')
        logger.error(f'Import error: {e}')
        return False

    except Exception as e:
        logger.error(f'Unexpected error during mesh visualization: {e}')
        logger.error(f'Error type: {type(e).__name__}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return False

def _handle_visualization_output(pPlotter, sFilename: Optional[str], iFlag_verbose: bool = False) -> bool:
    """
    Handle visualization output (save file or show interactive).

    Args:
        pPlotter: GeoVista plotter instance
        sFilename: Output filename or None for interactive
        iFlag_verbose: Enable verbose logging

    Returns:
        bool: True if successful, False otherwise
    """
    try:
        if sFilename is not None:
            # Save screenshot
            pPlotter.screenshot(sFilename)
            if iFlag_verbose:
                logger.info(f'✓ Visualization saved to: {sFilename}')

                # Verify file was created
                if os.path.exists(sFilename):
                    iFile_size = os.path.getsize(sFilename)
                    logger.info(f'  File size: {iFile_size / 1024:.1f} KB')
                else:
                    logger.warning(f'Screenshot command executed but file not found: {sFilename}')

            pPlotter.close()
            return True
        else:
            # Interactive display
            if iFlag_verbose:
                logger.info('Opening interactive visualization window...')
            pPlotter.show()
            return True

    except Exception as e:
        logger.error(f'Failed to handle visualization output: {e}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        try:
            pPlotter.close()
        except Exception:
            pass
        return False

def visualize_raster(self,
                    sFilename_out: Optional[str] = None,
                    dLongitude_focus_in: float = 0.0,
                    dLatitude_focus_in: float = 0.0,
                    dZoom_factor: float = 0.7,
                    iFlag_show_coastlines: bool = True,
                    iFlag_show_graticule: bool = True,
                    sColormap: str = 'viridis',
                    sCoastline_color: str = 'black',
                    dCoastline_width: float = 1.0,
                    iFlag_verbose: bool = False) -> bool:
    """
    Visualize source raster data using GeoVista.

    Creates 3D visualization of raster data by converting rasters to mesh format
    and displaying them with proper geographic context.

    Args:
        sFilename_out: Output screenshot file path. If None, displays interactive viewer.
        dLongitude_focus_in: Camera focal point longitude in degrees (-180 to 180).
        dLatitude_focus_in: Camera focal point latitude in degrees (-90 to 90).
        dZoom_factor: Camera zoom level. Higher values zoom in.
        iFlag_show_coastlines: Show coastline overlay.
        iFlag_show_graticule: Show coordinate grid with labels.
        sColormap: Matplotlib colormap name for raster visualization.
        sCoastline_color: Color for coastlines. Default is 'black'.
        dCoastline_width: Line width for coastlines. Default is 1.0.
        iFlag_verbose: If True, print detailed progress messages.

    Returns:
        True if visualization successful, False otherwise

    Note:
        - Requires 'geovista' package: pip install geovista
        - Converts raster data to mesh format for 3D visualization
        - Multiple rasters can be overlaid
    """
    if not self.aFilename_source_raster:
        logger.error('No source raster files available for visualization')
        return False

    if not _validate_output_path(sFilename_out):
        return False

    # Create configuration object
    config = VisualizationConfig(
        longitude_focus=dLongitude_focus_in,
        latitude_focus=dLatitude_focus_in,
        zoom_factor=dZoom_factor,
        show_coastlines=iFlag_show_coastlines,
        show_graticule=iFlag_show_graticule,
        colormap=sColormap,
        coastline_color=sCoastline_color,
        coastline_width=dCoastline_width,
        verbose=iFlag_verbose
    )

    try:
        import geovista as gv

        if config.verbose:
            logger.info(f'Visualizing {len(self.aFilename_source_raster)} raster file(s)...')

        # Setup plotter
        pPlotter = _setup_geovista_plotter(iFlag_off_screen=(sFilename_out is not None), iFlag_verbose=config.verbose)
        if pPlotter is None:
            return False

        # Process each raster file
        for idx, sFilename in enumerate(self.aFilename_source_raster, 1):
            if config.verbose:
                logger.info(f'Processing raster {idx}/{len(self.aFilename_source_raster)}: {os.path.basename(sFilename)}')

            try:
                pRaster = sraster(sFilename)
                pRaster.read_metadata()
                pRaster.create_raster_mesh()
                sFilename_raster_mesh = pRaster.sFilename_mesh

                if not os.path.exists(sFilename_raster_mesh):
                    logger.warning(f'Raster mesh file not found: {sFilename_raster_mesh}')
                    continue

                # Load raster mesh and add to visualization
                # This would require additional implementation to read the mesh
                # and extract raster values as cell data
                if config.verbose:
                    logger.info(f'  Created raster mesh: {sFilename_raster_mesh}')

                #we need to use the same apporoach as in visualize_source_mesh to load the mesh
                raster_mesh_info = rebuild_mesh_topology(sFilename_raster_mesh)
                if raster_mesh_info is None:
                    logger.warning(f'Failed to rebuild mesh topology for raster: {sFilename_raster_mesh}')
                    continue

                # Extract mesh data
                aVertex_longitude = raster_mesh_info['aVertex_longitude']
                aVertex_latitude = raster_mesh_info['aVertex_latitude']
                aConnectivity = raster_mesh_info['aConnectivity']
                aCellID = raster_mesh_info['aCellID']
                mesh = gv.Transform.from_unstructured(
                    aVertex_longitude,
                    aVertex_latitude,
                    connectivity=aConnectivity,
                    crs=CRS
                )
                name = f'Raster {idx} Cell ID'
                mesh.cell_data[name] = aCellID  # Placeholder for actual raster values
                sargs = {
                    "title": name,
                    "shadow": True,
                    "title_font_size": 10,
                    "label_font_size": 10,
                    "fmt": "%.0f",
                    "n_labels": 5,
                }
                pPlotter.add_mesh(mesh, scalar_bar_args=sargs)

            except Exception as e:
                logger.error(f"Error processing raster {sFilename}: {e}")
                continue

        # Configure camera and add geographic context
        _configure_camera(pPlotter, config)
        _add_geographic_context(pPlotter, config)

        # Handle output
        return _handle_visualization_output(pPlotter, sFilename_out, config.verbose)

    except ImportError as e:
        logger.error('GeoVista library not available. Install with: pip install geovista')
        logger.error(f'Import error: {e}')
        return False

    except Exception as e:
        logger.error(f'Unexpected error during raster visualization: {e}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return False

def _extract_target_mesh_data(sFilename: str, sVariable: str, iFlag_verbose: bool = False) -> Optional[Tuple[np.ndarray, int]]:
    """
    Extract variable data from target mesh file.

    Args:
        sFilename: Path to target mesh file
        sVariable: Variable name to extract
        iFlag_verbose: Enable verbose logging

    Returns:
        Tuple of (data_array, feature_count) or None if failed
    """
    try:
        # Open target mesh file
        pDataset = ogr.Open(sFilename, 0)  # Read-only
        if pDataset is None:
            logger.error(f'Failed to open target mesh file: {sFilename}')
            return None

        # Get first layer
        pLayer = pDataset.GetLayer(0)
        if pLayer is None:
            logger.error('Failed to get layer from target mesh dataset')
            pDataset = None
            return None

        # Get layer definition
        pLayerDefn = pLayer.GetLayerDefn()
        if pLayerDefn is None:
            logger.error('Failed to get layer definition from target mesh')
            pDataset = None
            return None

        # Get field information
        iFieldCount = pLayerDefn.GetFieldCount()
        nFeatures = pLayer.GetFeatureCount()

        if iFlag_verbose:
            logger.info(f'Target mesh contains {nFeatures} features with {iFieldCount} fields')

        # Check if variable field exists
        aField_names = [pLayerDefn.GetFieldDefn(i).GetName() for i in range(iFieldCount)]
        if sVariable not in aField_names:
            logger.error(f'Variable "{sVariable}" not found in target mesh')
            logger.error(f'Available fields: {", ".join(aField_names)}')
            pDataset = None
            return None

        if iFlag_verbose:
            logger.info(f'Extracting variable: {sVariable}')

        # Extract variable data from features, handling multipolygons correctly
        aData_list = []
        pLayer.ResetReading()
        pFeature = pLayer.GetNextFeature()
        iFeature_count = 0

        while pFeature is not None:
            pGeometry = pFeature.GetGeometryRef()
            if pGeometry is not None:
                sGeometry_type = pGeometry.GetGeometryName()

                if sGeometry_type == 'POLYGON':
                    # Single polygon - add one data value
                    try:
                        dField_value = pFeature.GetField(sVariable)
                        aData_list.append(dField_value if dField_value is not None else np.nan)
                    except Exception as e:
                        logger.warning(f'Error reading field {sVariable} from feature {iFeature_count}: {e}')
                        aData_list.append(np.nan)

                elif sGeometry_type == 'MULTIPOLYGON':
                    # Multipolygon - add the same data value for each polygon part
                    try:
                        dField_value = pFeature.GetField(sVariable)
                        dData_value = dField_value if dField_value is not None else np.nan

                        # Add the same data value for each polygon part in the multipolygon
                        nGeometryParts = pGeometry.GetGeometryCount()
                        iValid_parts = 0

                        for iPart in range(nGeometryParts):
                            pPolygon_part = pGeometry.GetGeometryRef(iPart)
                            if pPolygon_part is not None and pPolygon_part.IsValid():
                                aData_list.append(dData_value)
                                iValid_parts += 1

                        if iValid_parts == 0:
                            logger.warning(f'No valid parts found in multipolygon feature {iFeature_count}')
                            aData_list.append(np.nan)

                    except Exception as e:
                        logger.warning(f'Error reading field {sVariable} from multipolygon feature {iFeature_count}: {e}')
                        aData_list.append(np.nan)
                else:
                    logger.warning(f'Feature {iFeature_count} has unsupported geometry type: {sGeometry_type}')
                    aData_list.append(np.nan)
            else:
                logger.warning(f'Feature {iFeature_count} has no geometry')
                aData_list.append(np.nan)

            iFeature_count += 1
            pFeature = pLayer.GetNextFeature()

        # Close dataset
        pDataset = None

        if not aData_list:
            logger.error('No data extracted from target mesh')
            return None

        # Convert to numpy array
        aData = np.array(aData_list, dtype=np.float64)

        return aData, nFeatures

    except Exception as e:
        logger.error(f'Error extracting target mesh data: {e}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return None

def _create_target_mesh(pUraster_instance, aData: np.ndarray, sVariable: str, iFlag_verbose: bool = False) -> Optional[Tuple[Any, str, np.ndarray]]:
    """
    Create GeoVista mesh from uraster instance and attach data.

    Args:
        pUraster_instance: The uraster instance
        aData: Data array to attach to mesh
        sVariable: Variable name for the data
        iFlag_verbose: Enable verbose logging

    Returns:
        Tuple of (mesh, scalars_name, valid_cell_indices) or None if failed
    """
    try:
        import geovista as gv

        # Create masked connectivity array
        aConnectivity_masked = np.ma.masked_where(
            pUraster_instance.aConnectivity == -1,
            pUraster_instance.aConnectivity
        )

        # Validate connectivity indices
        aValid_connectivity = pUraster_instance.aConnectivity[pUraster_instance.aConnectivity >= 0]
        if len(aValid_connectivity) > 0:
            iMax_vertex_idx = len(pUraster_instance.aVertex_longititude) - 1
            if np.max(aValid_connectivity) > iMax_vertex_idx:
                logger.error(f'Connectivity contains invalid vertex index: '
                           f'max={np.max(aValid_connectivity)}, vertices={len(pUraster_instance.aVertex_longititude)}')
                return None

        # Transform to GeoVista unstructured mesh
        if iFlag_verbose:
            logger.info('Creating GeoVista mesh...')
        pMesh = gv.Transform.from_unstructured(
            pUraster_instance.aVertex_longititude,
            pUraster_instance.aVertex_latitude,
            connectivity=aConnectivity_masked,
            crs=CRS
        )

        if iFlag_verbose:
            logger.info(f'Created mesh with {pMesh.n_cells} cells and {pMesh.n_points} points')

        # Attach data to mesh
        sScalars = sVariable.capitalize()

        # Validate data array length matches mesh cells
        if len(aData) != pMesh.n_cells:
            logger.error(f'Data array length ({len(aData)}) does not match mesh cells ({pMesh.n_cells})')
            logger.error(f'This indicates a mismatch between mesh topology and extracted data')
            return None

        pMesh.cell_data[sScalars] = aData

        # Get valid cell indices (non-NaN values)
        aValid_data_mask = np.isfinite(aData)
        iN_valid = int(np.count_nonzero(aValid_data_mask))

        if iN_valid == 0:
            logger.warning(f'No valid cells to plot for variable "{sScalars}"')
            return None

        aValid_cell_indices = np.where(aValid_data_mask)[0]

        if iFlag_verbose:
            logger.info(f'Attached data "{sScalars}" to mesh cells')
            logger.info(f'Valid cells for visualization: {iN_valid}/{len(aData)}')

        return pMesh, sScalars, aValid_cell_indices

    except Exception as e:
        logger.error(f'Error creating target mesh: {e}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return None

def visualize_target_mesh(self,
                         sVariable_in: Optional[str] = None,
                         sUnit_in: Optional[str] = None,
                         sFilename_out: Optional[str] = None,
                         dLongitude_focus_in: Optional[float] = 0.0,
                         dLatitude_focus_in: Optional[float] = 0.0,
                         dZoom_factor: Optional[float] = 0.7,
                         iFlag_show_coastlines: Optional[bool] = True,
                         iFlag_show_graticule: Optional[bool] = True,
                         sColormap: Optional[str] = 'viridis',
                         sCoastline_color: Optional[str] = 'black',
                         dCoastline_width: Optional[float] = 1.0,
                         iFlag_create_animation: Optional[bool] = False,
                         iAnimation_frames: Optional[int] = 36,
                         dAnimation_speed: Optional[float] = 1.0,
                         sAnimation_format: Optional[str] = 'mp4',
                         iFlag_verbose: Optional[bool] = False) -> bool:
    """
    Visualize the target mesh with computed zonal statistics using GeoVista 3D rendering.

    Creates an interactive or saved 3D visualization of the mesh with cells colored
    by computed statistics (mean, min, max, std) from raster processing. Can also
    create rotating animations by generating multiple frames.

    Args:
        sVariable_in (str): Variable field name to visualize.
            Common values: 'mean', 'min', 'max', 'std', 'area'
        sUnit_in (str, optional): Unit label for the colorbar (e.g., 'mm', 'kg/m²').
            Default is empty string.
        sFilename_out (str, optional): Output screenshot file path.
            If None, displays interactive viewer. Supports: .png, .jpg, .svg
            For animations, this becomes the base filename (e.g., 'animation.mp4')
        dLongitude_focus_in (float, optional): Camera focal point longitude in degrees.
            Valid range: -180 to 180. Default is 0.0. For animations, this is the starting longitude.
        dLatitude_focus_in (float, optional): Camera focal point latitude in degrees.
            Valid range: -90 to 90. Default is 0.0.
        dZoom_factor (float, optional): Camera zoom level.
            Higher values zoom in. Default is 0.75.
        iFlag_show_coastlines (bool, optional): Show coastline overlay.
            Default is True.
        iFlag_show_graticule (bool, optional): Show coordinate grid with labels.
            Default is True.
        sColormap (str, optional): Matplotlib colormap name.
            Default is 'viridis'. Examples: 'plasma', 'coolwarm', 'jet', 'RdYlBu'
        sCoastline_color (str, optional): Color for coastlines. Default is 'black'.
            Examples: 'white', 'red', 'blue', 'gray', or RGB tuples like (1.0, 0.0, 0.0).
        dCoastline_width (float, optional): Line width for coastlines. Default is 1.0.
        iFlag_create_animation (bool, optional): Create rotating animation.
            Default is False. When True, generates frames for 360° rotation.
        iAnimation_frames (int, optional): Number of frames for 360° rotation.
            Default is 36 (10° per frame). More frames = smoother animation.
        dAnimation_speed (float, optional): Animation speed in degrees per frame.
            Default is 1.0. Calculated as 360 / iAnimation_frames if not specified.
        sAnimation_format (str, optional): Animation output format.
            Default is 'mp4'. Supports: 'mp4', 'gif', 'avi'
        iFlag_verbose (bool, optional): If True, print detailed progress messages.
            If False, only print error messages. Default is False.

    Returns:
        bool: True if visualization successful, False otherwise

    Raises:
        ImportError: If geovista package is not installed
        ValueError: If target mesh file or required data is not available

    Note:
        - Requires 'geovista' package: pip install geovista
        - Target mesh file must exist (created by run_remap method)
        - Specified variable must exist as a field in the target mesh
        - Interactive mode requires display environment
        - Animation mode requires 'imageio' package for video creation: pip install imageio[ffmpeg]
    """
    # Validate inputs
    if not self.sFilename_target_mesh:
        logger.error('No target mesh filename configured')
        return False

    if not os.path.exists(self.sFilename_target_mesh):
        logger.error(f'Target mesh file does not exist: {self.sFilename_target_mesh}')
        return False

    if not _validate_mesh_data(self):
        return False

    if not _validate_output_path(sFilename_out):
        return False

    # Set default variable if not provided
    sVariable = sVariable_in if sVariable_in else 'mean'
    if not isinstance(sVariable, str):
        logger.error(f'Variable name must be a string, got {type(sVariable).__name__}')
        return False

    # Create configuration objects
    config = VisualizationConfig(
        longitude_focus=dLongitude_focus_in,
        latitude_focus=dLatitude_focus_in,
        zoom_factor=dZoom_factor,
        show_coastlines=iFlag_show_coastlines,
        show_graticule=iFlag_show_graticule,
        colormap=sColormap,
        coastline_color=sCoastline_color,
        coastline_width=dCoastline_width,
        verbose=iFlag_verbose
    )

    animation_config = AnimationConfig(
        frames=iAnimation_frames,
        speed=dAnimation_speed,
        format=sAnimation_format
    ) if iFlag_create_animation else None

    try:
        import geovista as gv

        if config.verbose:
            logger.info(f'Loading target mesh data from: {self.sFilename_target_mesh}')

        # Extract data from target mesh
        pData_result = _extract_target_mesh_data(self.sFilename_target_mesh, sVariable, config.verbose)
        if pData_result is None:
            return False

        aData, nFeatures = pData_result
        self.nCell_target = nFeatures

        # Validate feature count matches source
        if self.nCell_source > 0 and self.nCell_target != self.nCell_source:
            logger.warning(f'Feature count mismatch: target={self.nCell_target}, source={self.nCell_source}')

        # Validate data
        aValid_data_mask = np.isfinite(aData)
        iValid_data_count = np.sum(aValid_data_mask)

        if iValid_data_count == 0:
            logger.error(f'All values for variable "{sVariable}" are invalid (NaN/Inf)')
            return False

        if iValid_data_count < len(aData):
            logger.warning(f'{len(aData) - iValid_data_count} of {len(aData)} values are invalid')

        # Log data statistics
        if config.verbose:
            aValid_values = aData[aValid_data_mask]
            logger.info(f'Data statistics for "{sVariable}":')
            logger.info(f'  - Valid values: {iValid_data_count}/{len(aData)}')
            logger.info(f'  - Min: {np.min(aValid_values):.4f}')
            logger.info(f'  - Max: {np.max(aValid_values):.4f}')
            logger.info(f'  - Mean: {np.mean(aValid_values):.4f}')
            logger.info(f'  - Std: {np.std(aValid_values):.4f}')

        # Create and validate mesh
        pMesh_result = _create_target_mesh(self, aData, sVariable, config.verbose)
        if pMesh_result is None:
            return False

        pMesh, sScalars, aValid_cell_indices = pMesh_result
        sUnit = sUnit_in if sUnit_in is not None else ""

        # Handle animation vs single frame visualization
        if animation_config is not None:
            return _handle_animation_visualization(
                pMesh, sScalars, aValid_cell_indices, sUnit, config, animation_config, sFilename_out
            )
        else:
            return _handle_single_frame_visualization(
                pMesh, sScalars, aValid_cell_indices, sUnit, config, sFilename_out
            )

    except ImportError as e:
        logger.error('GeoVista library not available. Install with: pip install geovista')
        logger.error(f'Import error: {e}')
        return False

    except Exception as e:
        logger.error(f'Unexpected error during target mesh visualization: {e}')
        logger.error(f'Error type: {type(e).__name__}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return False

def _handle_single_frame_visualization(pMesh, sScalars: str, aValid_cell_indices: np.ndarray,
                                     sUnit: str, pConfig: VisualizationConfig,
                                     sFilename: Optional[str]) -> bool:
    """
    Handle single frame visualization (static image or interactive).

    Args:
        pMesh: GeoVista mesh object
        sScalars: Name of scalar field to visualize
        aValid_cell_indices: Indices of valid cells to display
        sUnit: Unit string for colorbar
        pConfig: Visualization configuration
        sFilename: Output filename or None for interactive

    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Setup plotter
        pPlotter = _setup_geovista_plotter(iFlag_off_screen=(sFilename is not None), iFlag_verbose=pConfig.verbose)
        if pPlotter is None:
            return False

        # Configure scalar bar
        dSargs = {
            "title": f"{sScalars} / {sUnit}" if sUnit else sScalars,
            "shadow": True,
            "title_font_size": 12,
            "label_font_size": 10,
            "fmt": "%.2f",
            "n_labels": 5,
        }

        # Add mesh to plotter
        pMesh_valid = pMesh.extract_cells(aValid_cell_indices)
        pPlotter.add_mesh(pMesh_valid, scalars=sScalars, scalar_bar_args=dSargs, cmap=pConfig.colormap)

        # Configure camera and add geographic context
        _configure_camera(pPlotter, pConfig)
        _add_geographic_context(pPlotter, pConfig)

        # Handle output
        return _handle_visualization_output(pPlotter, sFilename, pConfig.verbose)

    except Exception as e:
        logger.error(f'Error in single frame visualization: {e}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return False

def _handle_animation_visualization(pMesh, sScalars: str, aValid_cell_indices: np.ndarray,
                                  sUnit: str, pConfig: VisualizationConfig,
                                  pAnimation_config: AnimationConfig, sFilename: str) -> bool:
    """
    Handle animation visualization.

    Args:
        pMesh: GeoVista mesh object
        sScalars: Name of scalar field to visualize
        aValid_cell_indices: Indices of valid cells to display
        sUnit: Unit string for colorbar
        pConfig: Visualization configuration
        pAnimation_config: Animation configuration
        sFilename: Output animation filename

    Returns:
        bool: True if successful, False otherwise
    """
    if sFilename is None:
        logger.error('Animation mode requires output filename')
        return False

    try:
        # Setup off-screen plotter for animation
        pPlotter = _setup_geovista_plotter(iFlag_off_screen=True, iFlag_verbose=pConfig.verbose)
        if pPlotter is None:
            return False

        # Configure scalar bar
        dSargs = {
            "title": f"{sScalars} / {sUnit}" if sUnit else sScalars,
            "shadow": True,
            "title_font_size": 12,
            "label_font_size": 10,
            "fmt": "%.2f",
            "n_labels": 5,
        }

        # Add mesh to plotter
        pMesh_valid = pMesh.extract_cells(aValid_cell_indices)
        pPlotter.add_mesh(pMesh_valid, scalars=sScalars, scalar_bar_args=dSargs, cmap=pConfig.colormap)

        # Configure initial camera position
        _configure_camera(pPlotter, pConfig)

        # Add geographic context
        _add_geographic_context(pPlotter, pConfig)

        # Reset the zoom factor to 1.0 so it won't zoom in too much during the animation
        pConfig.zoom_factor = 1.0

        # Create animation
        iFlag_success = _create_rotation_animation(
            pPlotter, sFilename, pConfig.longitude_focus, pConfig.latitude_focus,
            pConfig.zoom_factor, pAnimation_config.frames, pAnimation_config.speed,
            pAnimation_config.format, pConfig.verbose
        )

        pPlotter.close()
        return iFlag_success

    except Exception as e:
        logger.error(f'Error in animation visualization: {e}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return False

def _create_rotation_animation(pPlotter, sFilename_out, dLongitude_start, dLatitude_focus,
                               dZoom_factor, iAnimation_frames, dAnimation_speed, sAnimation_format, iFlag_verbose):
    """
    Create a rotating animation of the 3D globe visualization.

    Generates multiple frames by rotating the camera around the globe with enhanced
    camera movement patterns, then combines them into a video file.

    Args:
        pPlotter: GeoVista plotter instance with mesh already added
        sFilename_out (str): Output animation file path (e.g., 'animation.mp4')
        dLongitude_start (float): Starting longitude for rotation in degrees
        dLatitude_focus (float): Base latitude for camera focus in degrees
        dZoom_factor (float): Camera zoom level
        iAnimation_frames (int): Number of frames for 360° rotation
        dAnimation_speed (float): Degrees per frame
        sAnimation_format (str): Output format ('mp4', 'gif', 'avi')
        iFlag_verbose (bool): Enable verbose logging

    Returns:
        bool: True if animation created successfully, False otherwise
    """
    try:
        if iFlag_verbose:
            logger.info(f'Creating {iAnimation_frames}-frame rotation animation')
            logger.info(f'  - Starting longitude: {dLongitude_start:.1f}°')
            logger.info(f'  - Base latitude: {dLatitude_focus:.1f}°')
            logger.info(f'  - Rotation speed: {dAnimation_speed:.1f}°/frame')
            logger.info(f'  - Output format: {sAnimation_format}')

        # Animation parameters
        dEarth_radius = DEFAULT_EARTH_RADIUS
        dCamera_distance = dEarth_radius * DEFAULT_CAMERA_DISTANCE_MULTIPLIER
        dAmplitude_deg = 20.0  # Latitude oscillation amplitude
        dCycles = 1.0  # Number of sine cycles over full rotation
        dPhase = 0.0  # Phase shift for sine wave

        # Initialize movie recording
        pPlotter.open_movie(sFilename_out, framerate=30)

        if iFlag_verbose:
            logger.info('Generating animation frames...')

        for iFrame in range(iAnimation_frames):
            # Calculate current longitude with smooth rotation
            dLongitude_current = dLongitude_start + (iFrame * dAnimation_speed)
            dLongitude_current = dLongitude_current % 360.0  # Keep within [0, 360)
            if dLongitude_current > 180.0:
                dLongitude_current -= 360.0  # Convert to [-180, 180]

            # Enhanced latitude movement: sine-wave oscillation for dynamic viewing
            # This creates a more interesting camera path than fixed latitude
            dFrames_div = float(iAnimation_frames) if iAnimation_frames > 0 else 1.0
            dTheta = 2.0 * math.pi * (float(iFrame) / dFrames_div) * dCycles + dPhase
            dLatitude_current = float(dLatitude_focus) + dAmplitude_deg * math.sin(dTheta)

            # Clamp latitude to avoid pole singularities
            dLatitude_current = max(-89.9, min(89.9, dLatitude_current))

            # Convert to radians for calculations
            dLon_rad = math.radians(dLongitude_current)
            dLat_rad = math.radians(dLatitude_current)

            # Calculate focal point on Earth surface
            dX_focal = dEarth_radius * math.cos(dLat_rad) * math.cos(dLon_rad)
            dY_focal = dEarth_radius * math.cos(dLat_rad) * math.sin(dLon_rad)
            dZ_focal = dEarth_radius * math.sin(dLat_rad)

            # Calculate camera position away from Earth
            dX_camera = dCamera_distance * math.cos(dLat_rad) * math.cos(dLon_rad)
            dY_camera = dCamera_distance * math.cos(dLat_rad) * math.sin(dLon_rad)
            dZ_camera = dCamera_distance * math.sin(dLat_rad)

            aFocal_point = [dX_focal, dY_focal, dZ_focal]
            aCamera_position = [dX_camera, dY_camera, dZ_camera]

            # Update camera with smooth transitions
            pPlotter.camera.focal_point = aFocal_point
            pPlotter.camera.position = aCamera_position
            pPlotter.camera.up = [0, 0, 1]  # Maintain Z-up orientation

            # Apply zoom factor for consistent view
            pPlotter.camera.zoom(dZoom_factor)

            # Ensure axes remain visible throughout animation
            try:
                pPlotter.add_axes()
            except Exception:
                pass  # Axes may already exist

            # Render the current frame
            pPlotter.render()

            try:
                pPlotter.write_frame()

                if iFlag_verbose and (iFrame + 1) % max(1, iAnimation_frames // 10) == 0:
                    dProgress = ((iFrame + 1) / iAnimation_frames) * 100
                    logger.info(f'  Progress: {dProgress:.0f}% ({iFrame + 1}/{iAnimation_frames} frames)')

            except Exception as e:
                logger.error(f'Failed to render frame {iFrame + 1}: {e}')
                try:
                    pPlotter.close()
                except Exception:
                    pass
                return False

        # Close movie recording
        try:
            pPlotter.close()
        except Exception as e:
            logger.warning(f'Error closing plotter: {e}')

        # Validate output file creation
        if not os.path.exists(sFilename_out):
            logger.error('Animation file was not created')
            return False

        # Log success information
        iFile_size = os.path.getsize(sFilename_out)
        if iFlag_verbose:
            logger.info(f'✓ Animation created successfully: {sFilename_out}')
            logger.info(f'  File size: {iFile_size / (1024*1024):.2f} MB')
            logger.info(f'  Frames: {iAnimation_frames}')
            logger.info(f'  Format: {sAnimation_format.upper()}')
            logger.info(f'  Duration: ~{iAnimation_frames / 30:.1f} seconds at 30 FPS')

        return True

    except Exception as e:
        logger.error(f'Unexpected error during animation creation: {e}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        try:
            pPlotter.close()
        except Exception:
            pass
        return False

def rebuild_mesh_topology(sFilename_mesh_in, iFlag_verbose=False, sField_unique_id=None):
    """
    Rebuild mesh topology from source mesh file by extracting vertices,
    connectivity, and centroids for unstructured mesh processing.

    Args:
        sFilename_mesh_in (str): Path to the source mesh file (GeoJSON, Shapefile, etc.)
        iFlag_verbose (bool, optional): If True, print detailed progress messages.
            If False, only print error messages. Default is False.
        sField_unique_id (str, optional): Field name for unique cell IDs. If None, uses first field or feature index.

    Returns:
        dict: Comprehensive mesh topology information with keys:
            - 'vertices_longitude': np.ndarray of unique vertex longitudes
            - 'vertices_latitude': np.ndarray of unique vertex latitudes
            - 'connectivity': np.ndarray connectivity matrix
            - 'cell_centroids_longitude': np.ndarray of cell centroid longitudes
            - 'cell_centroids_latitude': np.ndarray of cell centroid latitudes
            - 'cell_ids': np.ndarray of cell IDs
            - 'area_min': float minimum cell area
            - 'area_max': float maximum cell area
            - 'area_mean': float mean cell area
            - 'max_vertices_per_cell': int maximum vertices per cell
            - 'num_cells': int total number of cells
            - 'num_vertices': int total number of unique vertices
            - 'success': bool whether processing was successful
        Returns None on failure.
    """
    try:
        # Open the input data source
        pDataset = ogr.Open(sFilename_mesh_in, 0)  # Read-only
        if pDataset is None:
            logger.error(f'Failed to open file: {sFilename_mesh_in}')
            return None
        if iFlag_verbose:
            logger.info(f'Successfully opened mesh file: {sFilename_mesh_in}')
        # Get the first layer
        pLayer = pDataset.GetLayer(0)
        if pLayer is None:
            logger.error('Failed to get layer from the dataset.')
            pDataset = None
            return None
        # Get layer information
        pLayerDefn = pLayer.GetLayerDefn()
        if pLayerDefn is None:
            logger.error('Failed to get layer definition.')
            pDataset = None
            return None
        nFeatures = pLayer.GetFeatureCount()
        nCell_source = nFeatures #if there is no invalid features, this is the number of cells
        iFieldCount = pLayerDefn.GetFieldCount()
        if nFeatures == 0:
            logger.warning('Layer contains no features.')
            pDataset = None
            return None
        aCellID= []  # Will be populated dynamically as features are processed
        if iFlag_verbose:
            logger.info(f'Processing {nFeatures} features with {iFieldCount} fields')
        # Get the first field name (assuming it contains the data variable)
        if sField_unique_id is None:
            sVariable = pLayerDefn.GetFieldDefn(0).GetName() if iFieldCount > 0 else None
            #search whether it has a field name used for id
            #for idx, pFeature_base in enumerate(pLayer):
            #    fid = pFeature_base.GetFID()
            #sVariable = None
        else:
            sVariable = sField_unique_id
        # Initialize lists for storing geometry data
        lons_list = []
        lats_list = []
        data_list = []
        area_list = []
        # Process features with enhanced error handling
        pLayer.ResetReading()
        iFeature_index = 0
        invalid_geometry_count = 0
        for pFeature in pLayer:
            if pFeature is None:
                continue
            pGeometry = pFeature.GetGeometryRef()
            if pGeometry is None:
                logger.warning(f'Feature {iFeature_index} has no geometry, skipping')
                iFeature_index += 1
                invalid_geometry_count += 1
                continue
            sGeometry_type = pGeometry.GetGeometryName()
            if sGeometry_type == 'POLYGON':
                try:
                    # Validate geometry before processing
                    if not pGeometry.IsValid():
                        logger.warning(f'Feature {iFeature_index} has invalid geometry')
                        invalid_geometry_count += 1
                        print(pGeometry.ExportToWkt())
                        continue
                    # Get coordinates of the polygon
                    aCoord = get_geometry_coordinates(pGeometry)
                    if aCoord is not None and len(aCoord) > 0:
                        # Validate coordinate bounds
                        lons = aCoord[:, 0]
                        lats = aCoord[:, 1]
                        # Check for reasonable coordinate ranges
                        if (np.any(lons < -180) or np.any(lons > 180) or
                            np.any(lats < -90) or np.any(lats > 90)):
                            logger.warning(f'Feature {iFeature_index} has coordinates outside valid range')
                        # Check for minimum polygon area (avoid degenerate polygons)
                        if len(aCoord) < 3:
                            logger.warning(f'Feature {iFeature_index} has fewer than 3 vertices, skipping')
                            iFeature_index += 1
                            invalid_geometry_count += 1
                            continue
                        lons_list.append(lons)
                        lats_list.append(lats)
                        # Calculate polygon area
                        try:
                            dArea = calculate_polygon_area(lons, lats)
                            area_list.append(dArea)
                        except Exception as area_error:
                            logger.warning(f'Could not calculate area for feature {iFeature_index}: {area_error}')
                            area_list.append(0.0)
                        # Get field data if available
                        if sVariable:
                            try:
                                field_value = pFeature.GetField(sVariable)
                                # Handle different field types
                                if field_value is not None:
                                    data_list.append(int(field_value))
                                else:
                                    data_list.append(0)
                                aCellID.append(int(field_value) if field_value is not None else iFeature_index)
                            except (ValueError, TypeError) as e:
                                logger.warning(f'Could not convert field value for feature {iFeature_index}: {e}')
                                data_list.append(iFeature_index)
                                aCellID.append(iFeature_index)
                        else:
                            data_list.append(iFeature_index)  # Use feature index as default
                            aCellID.append(iFeature_index)
                    else:
                        logger.warning(f'Failed to extract coordinates from feature {iFeature_index}')
                        invalid_geometry_count += 1
                except Exception as e:
                    logger.warning(f'Error processing feature {iFeature_index}: {str(e)}')
                    invalid_geometry_count += 1
            elif sGeometry_type == 'MULTIPOLYGON':
                try:
                    # Process multipolygon by extracting each constituent polygon
                    if iFlag_verbose:
                        logger.info(f'Processing multipolygon feature {iFeature_index} with {pGeometry.GetGeometryCount()} parts')
                    multipolygon_processed = False
                    for iPart in range(pGeometry.GetGeometryCount()):
                        pPolygon_part = pGeometry.GetGeometryRef(iPart)
                        if pPolygon_part is None:
                            logger.warning(f'Multipolygon part {iPart} is None in feature {iFeature_index}')
                            continue
                        if not pPolygon_part.IsValid():
                            logger.warning(f'Multipolygon part {iPart} has invalid geometry in feature {iFeature_index}')
                            continue
                        # Get coordinates of the polygon part
                        aCoord_part = get_geometry_coordinates(pPolygon_part)
                        if aCoord_part is not None and len(aCoord_part) > 0:
                            # Validate coordinate bounds for this part
                            lons_part = aCoord_part[:, 0]
                            lats_part = aCoord_part[:, 1]
                            # Check for reasonable coordinate ranges
                            if (np.any(lons_part < -180) or np.any(lons_part > 180) or
                                np.any(lats_part < -90) or np.any(lats_part > 90)):
                                logger.warning(f'Multipolygon part {iPart} in feature {iFeature_index} has coordinates outside valid range')
                            # Check for minimum polygon area (avoid degenerate polygons)
                            if len(aCoord_part) < 3:
                                logger.warning(f'Multipolygon part {iPart} in feature {iFeature_index} has fewer than 3 vertices, skipping part')
                                continue
                            lons_list.append(lons_part)
                            lats_list.append(lats_part)
                            # Calculate polygon area for this part
                            try:
                                dArea_part = calculate_polygon_area(lons_part, lats_part)
                                area_list.append(dArea_part)
                            except Exception as area_error:
                                logger.warning(f'Could not calculate area for multipolygon part {iPart} in feature {iFeature_index}: {area_error}')
                                area_list.append(0.0)
                            # For multipolygon, use the original feature index for all parts
                            # but track that this is a multipolygon part
                            if sVariable:
                                try:
                                    field_value = pFeature.GetField(sVariable)
                                    if field_value is not None:
                                        data_list.append(int(field_value))
                                    else:
                                        data_list.append(0)
                                    # Maintain original aCellID for multipolygon features
                                    # Each part gets the same CellID as the original feature
                                    aCellID.append(int(field_value) if field_value is not None else iFeature_index)
                                except (ValueError, TypeError) as e:
                                    logger.warning(f'Could not convert field value for multipolygon part {iPart} in feature {iFeature_index}: {e}')
                                    data_list.append(iFeature_index)
                                    aCellID.append(iFeature_index)
                            else:
                                data_list.append(iFeature_index)
                                aCellID.append(iFeature_index)
                            multipolygon_processed = True
                        else:
                            logger.warning(f'Failed to extract coordinates from multipolygon part {iPart} in feature {iFeature_index}')
                    if not multipolygon_processed:
                        logger.warning(f'No valid parts found in multipolygon feature {iFeature_index}')
                        invalid_geometry_count += 1
                    else:
                        if iFlag_verbose:
                            logger.info(f'Successfully processed multipolygon feature {iFeature_index}')
                except Exception as e:
                    logger.warning(f'Error processing multipolygon feature {iFeature_index}: {str(e)}')
                    invalid_geometry_count += 1
            elif sGeometry_type in ['POINT', 'LINESTRING']:
                logger.warning(f'Geometry type {sGeometry_type} not supported in feature {iFeature_index}, skipping')
                invalid_geometry_count += 1
            else:
                logger.warning(f'Unknown geometry type {sGeometry_type} in feature {iFeature_index}, skipping')
                invalid_geometry_count += 1
            iFeature_index += 1
        # Report processing statistics
        valid_mesh_cells = len(lons_list)
        if iFlag_verbose:
            logger.info(f'Feature processing summary:')
            logger.info(f'  - Total input features: {iFeature_index}')
            logger.info(f'  - Valid mesh cells created: {valid_mesh_cells}')
            logger.info(f'  - Invalid/skipped features: {invalid_geometry_count}')
            logger.info(f'  - Success rate: {((iFeature_index-invalid_geometry_count)/iFeature_index*100):.1f}%' if iFeature_index > 0 else '  - Success rate: 0%')
            # Report multipolygon handling statistics
            multipolygon_cells = valid_mesh_cells - (iFeature_index - invalid_geometry_count)
            if multipolygon_cells > 0:
                logger.info(f'  - Additional cells from multipolygons: {multipolygon_cells}')
                logger.info(f'  - Total mesh cells (including multipolygon parts): {valid_mesh_cells}')
        # Clean up dataset
        pDataset = None
        if not lons_list:
            logger.error('No valid polygon features found in mesh file')
            return None
        if iFlag_verbose:
            logger.info(f'Successfully processed {len(lons_list)} polygon features')
        # Calculate maximum vertices and pad coordinates efficiently
        try:
            if not lons_list:
                logger.error('No coordinate data found')
                return None
            max_vertices = max(len(coord) for coord in lons_list)
            if max_vertices == 0:
                logger.error('No vertices found in any polygon')
                return None
            nVertex_max = max_vertices
            if iFlag_verbose:
                logger.info(f'Maximum vertices per polygon: {max_vertices}')
            # Pre-allocate arrays for better memory efficiency
            num_polygons = len(lons_list)
            lons_padded = np.full((num_polygons, max_vertices), np.nan, dtype=np.float64)
            lats_padded = np.full((num_polygons, max_vertices), np.nan, dtype=np.float64)
            # Fill padded arrays efficiently
            for i, (lon_coords, lat_coords) in enumerate(zip(lons_list, lats_list)):
                # Ensure coordinates are numpy arrays with proper dtype
                lon_coords = np.asarray(lon_coords, dtype=np.float64)
                lat_coords = np.asarray(lat_coords, dtype=np.float64)
                # Validate coordinate data
                if len(lon_coords) != len(lat_coords):
                    logger.warning(f'Coordinate length mismatch in polygon {i}: lon={len(lon_coords)}, lat={len(lat_coords)}')
                    min_len = min(len(lon_coords), len(lat_coords))
                    lon_coords = lon_coords[:min_len]
                    lat_coords = lat_coords[:min_len]
                # Check for valid coordinate values
                if not (np.all(np.isfinite(lon_coords)) and np.all(np.isfinite(lat_coords))):
                    logger.warning(f'Invalid coordinates found in polygon {i}')
                    # Remove invalid coordinates
                    valid_mask = np.isfinite(lon_coords) & np.isfinite(lat_coords)
                    lon_coords = lon_coords[valid_mask]
                    lat_coords = lat_coords[valid_mask]
                coord_len = len(lon_coords)
                if coord_len > 0:
                    lons_padded[i, :coord_len] = lon_coords
                    lats_padded[i, :coord_len] = lat_coords
                else:
                    logger.warning(f'No valid coordinates remaining for polygon {i}')
            # Convert to the expected format for backward compatibility
            lons = lons_padded
            lats = lats_padded
        except Exception as e:
            logger.error(f'Error during coordinate padding: {str(e)}')
            logger.error(f'Traceback: {traceback.format_exc()}')
            return None
        # Calculate centroids efficiently using vectorized operations
        try:
            cell_lons_1d = []
            cell_lats_1d = []
            # Pre-allocate arrays for better performance
            cell_lons_1d = np.zeros(len(lons_list), dtype=np.float64)
            cell_lats_1d = np.zeros(len(lons_list), dtype=np.float64)
            for i in range(len(lons_list)):
                # Calculate centroid of each cell (ignoring NaN values)
                valid_mask = ~np.isnan(lons[i])
                if np.any(valid_mask):
                    valid_lons = lons[i][valid_mask]
                    valid_lats = lats[i][valid_mask]
                    # Use vectorized operations for better performance
                    centroid_lon = np.mean(valid_lons)
                    centroid_lat = np.mean(valid_lats)
                    # Validate centroid coordinates
                    if np.isfinite(centroid_lon) and np.isfinite(centroid_lat):
                        cell_lons_1d[i] = centroid_lon
                        cell_lats_1d[i] = centroid_lat
                    else:
                        logger.warning(f'Invalid centroid calculated for cell {i}: lon={centroid_lon}, lat={centroid_lat}')
                        # Use geometric center of bounding box as fallback
                        if len(valid_lons) > 0 and len(valid_lats) > 0:
                            cell_lons_1d[i] = (np.min(valid_lons) + np.max(valid_lons)) / 2.0
                            cell_lats_1d[i] = (np.min(valid_lats) + np.max(valid_lats)) / 2.0
                        else:
                            cell_lons_1d[i] = 0.0
                            cell_lats_1d[i] = 0.0
                else:
                    logger.warning(f'No valid coordinates found for cell {i}')
                    cell_lons_1d[i] = 0.0
                    cell_lats_1d[i] = 0.0
            if iFlag_verbose:
                logger.info(f'Calculated centroids for {len(cell_lons_1d)} cells')
            # Validate centroid ranges
            lon_range = (np.min(cell_lons_1d), np.max(cell_lons_1d))
            lat_range = (np.min(cell_lats_1d), np.max(cell_lats_1d))
            if not (-180 <= lon_range[0] <= 180 and -180 <= lon_range[1] <= 180):
                logger.warning(f'Longitude centroids outside valid range: {lon_range}')
            if not (-90 <= lat_range[0] <= 90 and -90 <= lat_range[1] <= 90):
                logger.warning(f'Latitude centroids outside valid range: {lat_range}')
        except Exception as e:
            logger.error(f'Error during centroid calculation: {str(e)}')
            return None
        # Extract unique vertices and connectivity
        try:
            if iFlag_verbose:
                logger.info('Extracting unique vertices and connectivity...')
            xv, yv, connectivity, vertex_to_index = extract_unique_vertices_and_connectivity(
                lons_list, lats_list
            )
            if xv is None or yv is None or connectivity is None:
                logger.error('Failed to extract unique vertices and connectivity')
                return None
            if iFlag_verbose:
                logger.info(f'Extracted {len(xv)} unique vertices')
                logger.info(f'Created connectivity matrix with shape: {connectivity.shape}')
        except Exception as e:
            logger.error(f'Error during vertex/connectivity extraction: {str(e)}')
            return None
        # Store results in class attributes
        aVertex_longititude = xv
        aVertex_latitude = yv
        aCenter_longititude = cell_lons_1d
        aCenter_latitude = cell_lats_1d
        aConnectivity = connectivity
        # Ensure aCellID matches the number of valid mesh cells
        if len(aCellID) != len(cell_lons_1d):
            logger.warning(f"aCellID length ({len(aCellID)}) doesn't match mesh cells ({len(cell_lons_1d)})")
            if len(aCellID) > len(cell_lons_1d):
                # Truncate aCellID to match mesh cells
                logger.warning("Truncating aCellID to match mesh cell count")
                aCellID = aCellID[:len(cell_lons_1d)]
            else:
                # Extend aCellID with sequential indices
                logger.warning("Extending aCellID with sequential indices to match mesh cell count")
                missing_count = len(cell_lons_1d) - len(aCellID)
                aCellID.extend(range(len(aCellID), len(aCellID) + missing_count))
        aCellID = np.array(aCellID)
        if iFlag_verbose:
            logger.info(f'Final aCellID array length: {len(aCellID)}')
            logger.info(f'aCellID range: [{np.min(aCellID)}, {np.max(aCellID)}]')
        # Calculate and store area statistics
        if area_list:
            area_array = np.array(area_list)
            valid_areas = area_array[area_array > 0]  # Exclude zero areas from statistics
            if len(valid_areas) > 0:
                dArea_min = float(np.min(valid_areas))
                dArea_max = float(np.max(valid_areas))
                dArea_mean = float(np.mean(valid_areas))
                dArea_max = float(np.max(valid_areas))
                dArea_min = float(np.min(valid_areas))
                if iFlag_verbose:
                    logger.info(f'Mesh area statistics:')
                    logger.info(f'  - Min area: {dArea_min:.6f}')
                    logger.info(f'  - Max area: {dArea_max:.6f}')
                    logger.info(f'  - Mean area: {dArea_mean:.6f}')
            else:
                logger.warning('No valid polygon areas calculated')
                dArea_min = 0.0
                dArea_max = 0.0
                dArea_mean = 0.0
        # Enhanced validation of final results
        validation_passed = True
        if len(aVertex_longititude) == 0:
            logger.error('No unique vertices extracted')
            validation_passed = False
        if len(aCenter_longititude) != len(lons_list):
            logger.error(f'Centroid count mismatch: expected {len(lons_list)}, got {len(aCenter_longititude)}')
            validation_passed = False
        if aConnectivity is None or aConnectivity.size == 0:
            logger.error('Empty connectivity matrix')
            validation_passed = False
        # Validate connectivity indices
        if aConnectivity is not None:
            max_vertex_index = len(aVertex_longititude) - 1
            valid_connectivity = aConnectivity[aConnectivity >= 0]
            if len(valid_connectivity) > 0 and np.max(valid_connectivity) > max_vertex_index:
                logger.error('Connectivity matrix contains invalid vertex indices')
                validation_passed = False
        # Check for reasonable mesh bounds
        if len(aVertex_longititude) > 0:
            vertex_lon_range = (np.min(aVertex_longititude), np.max(aVertex_longititude))
            vertex_lat_range = (np.min(aVertex_latitude), np.max(aVertex_latitude))
            if not (-180 <= vertex_lon_range[0] <= 180 and -180 <= vertex_lon_range[1] <= 180):
                logger.warning(f'Vertex longitudes outside valid range: {vertex_lon_range}')
            if not (-90 <= vertex_lat_range[0] <= 90 and -90 <= vertex_lat_range[1] <= 90):
                logger.warning(f'Vertex latitudes outside valid range: {vertex_lat_range}')
        if not validation_passed:
            logger.error('Mesh topology rebuild failed validation')
            return None
        if iFlag_verbose:
            logger.info('Mesh topology successfully rebuilt')
            logger.info(f'Final mesh statistics:')
            logger.info(f'  - Unique vertices: {len(aVertex_longititude)}')
            logger.info(f'  - Mesh cells: {len(aCenter_longititude)}')
            logger.info(f'  - Max vertices per cell: {nVertex_max}')
            logger.info(f'  - Connectivity shape: {aConnectivity.shape}')
            logger.info(f'  - Vertex longitude range: [{np.min(aVertex_longititude):.3f}, {np.max(aVertex_longititude):.3f}]')
            logger.info(f'  - Vertex latitude range: [{np.min(aVertex_latitude):.3f}, {np.max(aVertex_latitude):.3f}]')

        # Return comprehensive mesh topology information
        mesh_info = {
            'vertices_longitude': aVertex_longititude,
            'vertices_latitude': aVertex_latitude,
            'connectivity': aConnectivity,
            'cell_centroids_longitude': aCenter_longititude,
            'cell_centroids_latitude': aCenter_latitude,
            'cell_ids': aCellID,
            'area_min': dArea_min,
            'area_max': dArea_max,
            'area_mean': dArea_mean,
            'max_vertices_per_cell': nVertex_max,
            'num_cells': len(aCenter_longititude),
            'num_vertices': len(aVertex_longititude),
            'success': True
        }

        return mesh_info
    except Exception as e:
        logger.error(f'Unexpected error in rebuild_mesh_topology: {str(e)}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return None