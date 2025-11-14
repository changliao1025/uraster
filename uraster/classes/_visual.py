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
from uraster.classes.sraster import sraster
from uraster import utility
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
                raster_mesh_info = utility.rebuild_mesh_topology(
                    sFilename_raster_mesh,
                    sField_unique_id=self.sField_unique_id
                )
                if raster_mesh_info is None:
                    logger.warning(f'Failed to rebuild mesh topology for raster: {sFilename_raster_mesh}')
                    continue

                # Extract mesh data
                aVertex_longitude = raster_mesh_info['vertices_longitude']
                aVertex_latitude = raster_mesh_info['vertices_latitude']
                aConnectivity = raster_mesh_info['connectivity']
                aCellID = raster_mesh_info['cell_ids']
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

