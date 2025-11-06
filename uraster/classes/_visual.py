"""
Visualization module for uraster class.

This module contains all visualization-related methods that were moved from the main uraster class
to reduce the size of the main uraster.py file and improve code organization.
"""

import os
import logging
import traceback
import math
import numpy as np
from osgeo import gdal, ogr
from uraster.classes.sraster import sraster
gdal.UseExceptions()
# Set up logging
logger = logging.getLogger(__name__)
crs = "EPSG:4326"

def visualize_source_mesh(self,
                          sFilename_out=None,
                          dLongitude_focus_in=0.0,
                          dLatitude_focus_in=0.0,
                          dZoom_factor=0.7,
                          iFlag_show_coastlines=True,
                          iFlag_show_graticule=True,
                          iFlag_verbose=False):
    """
    Visualize the source mesh topology using GeoVista 3D globe rendering.

    Creates an interactive or saved 3D visualization of the unstructured mesh
    with proper geographic context including coastlines and coordinate grid.

    Args:
        sFilename_out (str, optional): Output screenshot file path.
            If None, displays interactive viewer. Supports formats: .png, .jpg, .svg
        dLongitude_focus_in (float, optional): Camera focal point longitude in degrees.
            Valid range: -180 to 180. Default is 0.0 (prime meridian).
        dLatitude_focus_in (float, optional): Camera focal point latitude in degrees.
            Valid range: -90 to 90. Default is 0.0 (equator).
        dZoom_factor (float, optional): Camera zoom level.
            Higher values zoom in. Default is 0.7.
        iFlag_show_coastlines (bool, optional): Show coastline overlay.
            Default is True.
        iFlag_show_graticule (bool, optional): Show coordinate grid with labels.
            Default is True.
        iFlag_verbose (bool, optional): If True, print detailed progress messages.
            If False, only print error messages. Default is False.

    Returns:
        bool: True if visualization successful, False otherwise

    Note:
        - Requires 'geovista' package: pip install geovista
        - Interactive mode requires display environment
        - Mesh topology must be built before visualization (call rebuild_mesh_topology first)
    """
    # Validate mesh data availability
    if self.aVertex_longititude is None or self.aVertex_latitude is None:
        logger.error('Mesh vertices not available. Build mesh topology first.')
        return False

    if self.aConnectivity is None:
        logger.error('Mesh connectivity not available. Build mesh topology first.')
        return False

    if len(self.aVertex_longititude) == 0 or len(self.aVertex_latitude) == 0:
        logger.error('Mesh vertices are empty.')
        return False

    if self.aConnectivity.size == 0:
        logger.error('Mesh connectivity is empty.')
        return False

    # Validate focus coordinates
    dLongitude_focus = dLongitude_focus_in if dLongitude_focus_in is not None else 0.0
    dLatitude_focus = dLatitude_focus_in if dLatitude_focus_in is not None else 0.0

    if not (-180 <= dLongitude_focus <= 180):
        logger.warning(f'Longitude focus {dLongitude_focus} out of range [-180, 180], clamping')
        dLongitude_focus = np.clip(dLongitude_focus, -180, 180)

    if not (-90 <= dLatitude_focus <= 90):
        logger.warning(f'Latitude focus {dLatitude_focus} out of range [-90, 90], clamping')
        dLatitude_focus = np.clip(dLatitude_focus, -90, 90)

    # Validate zoom factor
    if dZoom_factor <= 0:
        logger.warning(f'Invalid zoom factor {dZoom_factor}, using default 0.7')
        dZoom_factor = 0.7

    # Validate output file path if provided
    if sFilename_out is not None:
        if not isinstance(sFilename_out, str) or not sFilename_out.strip():
            logger.error('Output filename must be a non-empty string')
            return False

        # Check output directory exists
        output_dir = os.path.dirname(sFilename_out)
        if output_dir and not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir, exist_ok=True)
                if iFlag_verbose:
                    logger.info(f'Created output directory: {output_dir}')
            except Exception as e:
                logger.error(f'Cannot create output directory {output_dir}: {e}')
                return False

        # Check supported file extensions
        valid_extensions = ['.png', '.jpg', '.jpeg', '.svg', '.tif', '.tiff']
        file_ext = os.path.splitext(sFilename_out)[1].lower()
        if file_ext not in valid_extensions:
            logger.warning(f'File extension {file_ext} may not be supported. Recommended: .png, .jpg, .svg')

    # Import geovista with error handling
    try:
        import geovista as gv
        if iFlag_verbose:
            logger.info('GeoVista library imported successfully')
    except ImportError as e:
        logger.error('GeoVista library not available. Install with: pip install geovista')
        logger.error(f'Import error: {e}')
        return False

    try:
        if iFlag_verbose:
            logger.info('Creating mesh visualization...')
            logger.info(f'  - Vertices: {len(self.aVertex_longititude)}')
            logger.info(f'  - Connectivity shape: {self.aConnectivity.shape}')
            logger.info(f'  - Focus: ({dLongitude_focus:.2f}°, {dLatitude_focus:.2f}°)')
            logger.info(f'  - Zoom factor: {dZoom_factor}')

        # Prepare mesh metadata
        name = 'Mesh Cell ID'
        sUnit = ""

        # Validate connectivity array structure
        if self.aConnectivity.size == 0:
            logger.error('Connectivity array is empty')
            return False

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
                logger.error(f'Connectivity contains invalid vertex index: max={np.max(valid_connectivity)}, vertices={len(self.aVertex_longititude)}')
                return False

        # Transform to GeoVista unstructured mesh
        mesh = gv.Transform.from_unstructured(
            self.aVertex_longititude,
            self.aVertex_latitude,
            connectivity=connectivity_masked,
            crs=crs
        )

        # Validate cell data array length matches mesh cells
        if len(self.aCellID) != mesh.n_cells:
            logger.error(f'Cell ID array length ({len(self.aCellID)}) does not match mesh cells ({mesh.n_cells})')
            logger.error(f'This indicates a mismatch between mesh topology and cell data')
            return False

        mesh.cell_data[name] = self.aCellID

        if iFlag_verbose:
            logger.info(f'Created GeoVista mesh with {mesh.n_cells} cells and {mesh.n_points} points')

        # Create 3D plotter with enhanced error handling
        try:
            if sFilename_out is not None:
                # Off-screen rendering for saving files
                plotter = gv.GeoPlotter(off_screen=True)
                if iFlag_verbose:
                    logger.debug('Created off-screen plotter for file output')
            else:
                # Interactive rendering
                plotter = gv.GeoPlotter()
                if iFlag_verbose:
                    logger.debug('Created interactive plotter')

        except Exception as e:
            logger.error(f'Failed to create GeoVista plotter: {e}')
            logger.error('This may be due to missing graphics context or display')
            if sFilename_out is not None:
                logger.error('For headless systems, ensure proper OpenGL/Mesa setup')
            else:
                logger.error('For interactive mode, ensure display environment (X11/Wayland) is available')
            return False

        # Configure scalar bar (colorbar) appearance
        sargs = {
            "title": f"{name} / {sUnit}" if sUnit else name,
            "shadow": True,
            "title_font_size": 10,
            "label_font_size": 10,
            "fmt": "%.0f",  # Integer formatting for cell IDs
            "n_labels": 5,
        }

        # Add mesh to plotter
        plotter.add_mesh(mesh, scalars=name, scalar_bar_args=sargs)
        # Configure camera position and focus
        try:
            # Use PyVista/VTK coordinate conversion instead of deprecated geodesic


            # Convert longitude/latitude to radians
            lon_rad = math.radians(dLongitude_focus)
            lat_rad = math.radians(dLatitude_focus)

            # Earth radius (approximately 6371 km, but use normalized units)
            earth_radius = 1.0
            camera_distance = earth_radius * 3.0  # Position camera 3x earth radius away

            # Convert spherical coordinates to Cartesian (x, y, z)
            # Focal point is ON the Earth surface at the specified coordinates
            x_focal = earth_radius * math.cos(lat_rad) * math.cos(lon_rad)
            y_focal = earth_radius * math.cos(lat_rad) * math.sin(lon_rad)
            z_focal = earth_radius * math.sin(lat_rad)

            # Camera position is AWAY from Earth at the same angular position
            # This creates a proper viewing angle from outside looking in
            x_camera = camera_distance * math.cos(lat_rad) * math.cos(lon_rad)
            y_camera = camera_distance * math.cos(lat_rad) * math.sin(lon_rad)
            z_camera = camera_distance * math.sin(lat_rad)

            focal_point = [x_focal, y_focal, z_focal]
            camera_position = [x_camera, y_camera, z_camera]

            # Validate camera setup - ensure camera and focal point are different
            camera_distance_check = math.sqrt(
                (x_camera - x_focal)**2 + (y_camera - y_focal)**2 + (z_camera - z_focal)**2
            )

            if camera_distance_check < 0.1:  # Too close or identical
                logger.warning('Camera and focal point are too close, using default view')
                raise ValueError('Invalid camera positioning')

            plotter.camera.focal_point = focal_point
            plotter.camera.position = camera_position
            plotter.camera.zoom(dZoom_factor)

            # Set up vector for proper orientation
            plotter.camera.up = [0, 0, 1]  # Z-up orientation

            # Verify camera setup was successful
            if hasattr(plotter.camera, 'position') and plotter.camera.position:
                if iFlag_verbose:
                    logger.debug(f'Camera configured successfully:')
                    logger.debug(f'  Focal point: [{x_focal:.3f}, {y_focal:.3f}, {z_focal:.3f}]')
                    logger.debug(f'  Camera position: [{x_camera:.3f}, {y_camera:.3f}, {z_camera:.3f}]')
                    logger.debug(f'  Distance: {camera_distance_check:.3f}')
            else:
                logger.warning('Camera configuration may not have been applied correctly')

        except Exception as e:
            logger.warning(f'Error setting camera position: {e}. Using default view.')
            # Reset to default camera view
            try:
                plotter.reset_camera()
            except Exception:
                pass

        # Add geographic context
        if iFlag_show_coastlines:
            try:
                plotter.add_coastlines()
                if iFlag_verbose:
                    logger.debug('Added coastlines overlay')
            except Exception as e:
                logger.warning(f'Could not add coastlines: {e}')

        # Add coordinate axes
        try:
            plotter.add_axes()
            if iFlag_verbose:
                logger.debug('Added coordinate axes')
        except Exception as e:
            logger.warning(f'Could not add axes: {e}')

        # Add graticule (coordinate grid)
        if iFlag_show_graticule:
            try:
                plotter.add_graticule(show_labels=True)
                if iFlag_verbose:
                    logger.debug('Added coordinate graticule with labels')
            except Exception as e:
                logger.warning(f'Could not add graticule: {e}')

        # Output or display
        if sFilename_out is not None:
            # Save screenshot
            try:
                plotter.screenshot(sFilename_out)
                if iFlag_verbose:
                    logger.info(f'✓ Visualization saved to: {sFilename_out}')

                    # Verify file was created
                    if os.path.exists(sFilename_out):
                        file_size = os.path.getsize(sFilename_out)
                        logger.info(f'  File size: {file_size / 1024:.1f} KB')
                    else:
                        logger.warning(f'Screenshot command executed but file not found: {sFilename_out}')

                plotter.close()
                return True

            except Exception as e:
                logger.error(f'Failed to save screenshot: {e}')
                logger.error(f'Traceback: {traceback.format_exc()}')
                plotter.close()
                return False
        else:
            # Interactive display
            try:
                if iFlag_verbose:
                    logger.info('Opening interactive visualization window...')
                plotter.show()
                return True
            except Exception as e:
                logger.error(f'Failed to display interactive visualization: {e}')
                logger.error(f'Ensure display environment is available (X11, Wayland, etc.)')
                logger.error(f'Traceback: {traceback.format_exc()}')
                plotter.close()
                return False

    except ImportError as e:
        logger.error(f'Missing required GeoVista dependencies: {e}')
        return False

    except Exception as e:
        logger.error(f'Unexpected error during mesh visualization: {e}')
        logger.error(f'Error type: {type(e).__name__}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return False


def visualize_raster(self):
    """
    Visualize source raster data using GeoVista.

    Note:
        Not yet implemented. Placeholder for future raster visualization.
    """
    #it is possible to visualize the raster using geovista as well
    #if more than one rasters is used, then we may overlap them one by one
    for idx, sFilename in enumerate(self.aFilename_source_raster, 1):
        print(f"\n[{idx}] {sFilename}")
        try:
            pRaster = sraster(sFilename)
            pRaster.create_raster_mesh()
            sFilename_raster_mesh = pRaster.sFilename_mesh

            #use this mesh to visualize, with the raster value as the cell data

        except Exception as e:
            logger.error(f"Error reading raster info: {e}")
    return


def visualize_target_mesh(self, sVariable_in=None,
                           sUnit_in=None,
                           sFilename_out=None,
                           dLongitude_focus_in=0.0,
                           dLatitude_focus_in=0.0,
                           dZoom_factor=0.7,
                           iFlag_show_coastlines=True,
                           iFlag_show_graticule=True,
                           sColormap='viridis',
                           iFlag_create_animation=False,
                           iAnimation_frames=36,
                           dAnimation_speed=1.0,
                           sAnimation_format='mp4',
                           iFlag_verbose=False):
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
    # Validate target mesh file
    if not self.sFilename_target_mesh:
        logger.error('No target mesh filename configured')
        return False

    if not os.path.exists(self.sFilename_target_mesh):
        logger.error(f'Target mesh file does not exist: {self.sFilename_target_mesh}')
        return False

    # Validate mesh topology data
    if self.aVertex_longititude is None or self.aVertex_latitude is None:
        logger.error('Mesh vertices not available. Build mesh topology first.')
        return False

    if self.aConnectivity is None:
        logger.error('Mesh connectivity not available. Build mesh topology first.')
        return False

    if len(self.aVertex_longititude) == 0 or len(self.aVertex_latitude) == 0:
        logger.error('Mesh vertices are empty.')
        return False

    # Validate variable name
    if not sVariable_in:
        logger.warning('No variable specified, defaulting to "mean"')
        sVariable = 'mean'
    else:
        sVariable = sVariable_in

    if not isinstance(sVariable, str):
        logger.error(f'Variable name must be a string, got {type(sVariable).__name__}')
        return False

    # Validate focus coordinates
    dLongitude_focus = dLongitude_focus_in if dLongitude_focus_in is not None else 0.0
    dLatitude_focus = dLatitude_focus_in if dLatitude_focus_in is not None else 0.0

    if not (-180 <= dLongitude_focus <= 180):
        logger.warning(f'Longitude focus {dLongitude_focus} out of range [-180, 180], clamping')
        dLongitude_focus = np.clip(dLongitude_focus, -180, 180)

    if not (-90 <= dLatitude_focus <= 90):
        logger.warning(f'Latitude focus {dLatitude_focus} out of range [-90, 90], clamping')
        dLatitude_focus = np.clip(dLatitude_focus, -90, 90)

    # Validate zoom factor
    if dZoom_factor <= 0:
        logger.warning(f'Invalid zoom factor {dZoom_factor}, using default 1.0')
        dZoom_factor = 0.7

    # Validate output file path if provided
    if sFilename_out is not None:
        if not isinstance(sFilename_out, str) or not sFilename_out.strip():
            logger.error('Output filename must be a non-empty string')
            return False

        # Check output directory exists
        output_dir = os.path.dirname(sFilename_out)
        if output_dir and not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir, exist_ok=True)
                logger.info(f'Created output directory: {output_dir}')
            except Exception as e:
                logger.error(f'Cannot create output directory {output_dir}: {e}')
                return False

        # Check supported file extensions
        valid_extensions = ['.png', '.jpg', '.jpeg', '.svg', '.tif', '.tiff', '.mp4', '.gif', '.avi']
        file_ext = os.path.splitext(sFilename_out)[1].lower()
        if file_ext not in valid_extensions:
            logger.warning(f'File extension {file_ext} may not be supported. Recommended: .png, .jpg, .svg')

    # Import geovista with error handling
    try:
        import geovista as gv
        if iFlag_verbose:
            logger.info('GeoVista library imported successfully')
    except ImportError as e:
        logger.error('GeoVista library not available. Install with: pip install geovista')
        logger.error(f'Import error: {e}')
        return False

    try:
        if iFlag_verbose:
            logger.info(f'Loading target mesh data from: {self.sFilename_target_mesh}')

        # Open target mesh file
        pDataset = ogr.Open(self.sFilename_target_mesh, 0)  # Read-only
        if pDataset is None:
            logger.error(f'Failed to open target mesh file: {self.sFilename_target_mesh}')
            return False

        # Get first layer
        pLayer = pDataset.GetLayer(0)
        if pLayer is None:
            logger.error('Failed to get layer from target mesh dataset')
            pDataset = None
            return False

        # Get layer definition
        pLayerDefn = pLayer.GetLayerDefn()
        if pLayerDefn is None:
            logger.error('Failed to get layer definition from target mesh')
            pDataset = None
            return False

        # Get field information
        iFieldCount = pLayerDefn.GetFieldCount()
        nFeatures = pLayer.GetFeatureCount()
        self.nCell_target = nFeatures

        if iFlag_verbose:
            logger.info(f'Target mesh contains {nFeatures} features with {iFieldCount} fields')

        # Check if variable field exists
        field_names = [pLayerDefn.GetFieldDefn(i).GetName() for i in range(iFieldCount)]
        if sVariable not in field_names:
            logger.error(f'Variable "{sVariable}" not found in target mesh')
            logger.error(f'Available fields: {", ".join(field_names)}')
            pDataset = None
            return False

        if iFlag_verbose:
            logger.info(f'Extracting variable: {sVariable}')

        # Validate feature count matches source
        if self.nCell_source > 0 and self.nCell_target != self.nCell_source:
            logger.warning(f'Feature count mismatch: target={self.nCell_target}, source={self.nCell_source}')

        # Extract variable data from features, handling multipolygons correctly
        # This must match the logic used in rebuild_mesh_topology() where multipolygon
        # parts are added as separate mesh cells
        data_list = []
        pLayer.ResetReading()
        pFeature = pLayer.GetNextFeature()
        feature_count = 0

        while pFeature is not None:
            pGeometry = pFeature.GetGeometryRef()
            if pGeometry is not None:
                sGeometry_type = pGeometry.GetGeometryName()

                if sGeometry_type == 'POLYGON':
                    # Single polygon - add one data value
                    try:
                        field_value = pFeature.GetField(sVariable)
                        if field_value is not None:
                            data_list.append(field_value)
                        else:
                            # Handle None values as NaN
                            data_list.append(np.nan)
                            #logger.warning(f'Feature {feature_count} has None value for {sVariable}')
                    except Exception as e:
                        logger.warning(f'Error reading field {sVariable} from feature {feature_count}: {e}')
                        data_list.append(np.nan)

                elif sGeometry_type == 'MULTIPOLYGON':
                    # Multipolygon - add the same data value for each polygon part
                    # This matches the mesh topology building logic
                    try:
                        field_value = pFeature.GetField(sVariable)
                        data_value = field_value if field_value is not None else np.nan

                        if field_value is None:
                            #logger.warning(f'Feature {feature_count} has None value for {sVariable}')
                            pass

                        # Add the same data value for each polygon part in the multipolygon
                        nGeometryParts = pGeometry.GetGeometryCount()
                        valid_parts = 0

                        for iPart in range(nGeometryParts):
                            pPolygon_part = pGeometry.GetGeometryRef(iPart)
                            if pPolygon_part is not None and pPolygon_part.IsValid():
                                # Only add data for valid polygon parts (matching mesh topology logic)
                                data_list.append(data_value)
                                valid_parts += 1
                            else:
                                logger.warning(f'Invalid polygon part {iPart} in multipolygon feature {feature_count}')

                        if valid_parts == 0:
                            logger.warning(f'No valid parts found in multipolygon feature {feature_count}')
                            data_list.append(np.nan)  # Add at least one value to maintain count
                        else:
                            logger.debug(f'Added {valid_parts} data values for multipolygon feature {feature_count}')

                    except Exception as e:
                        logger.warning(f'Error reading field {sVariable} from multipolygon feature {feature_count}: {e}')
                        data_list.append(np.nan)
                else:
                    logger.warning(f'Feature {feature_count} has unsupported geometry type: {sGeometry_type}')
                    data_list.append(np.nan)
            else:
                logger.warning(f'Feature {feature_count} has no geometry')
                data_list.append(np.nan)

            feature_count += 1
            pFeature = pLayer.GetNextFeature()

        # Close dataset
        pDataset = None

        if not data_list:
            logger.error('No data extracted from target mesh')
            return False

        # Convert to numpy array
        data = np.array(data_list, dtype=np.float64)

        # Validate data
        valid_data_mask = np.isfinite(data)
        valid_data_count = np.sum(valid_data_mask)

        if valid_data_count == 0:
            logger.error(f'All values for variable "{sVariable}" are invalid (NaN/Inf)')
            return False

        if valid_data_count < len(data):
            logger.warning(f'{len(data) - valid_data_count} of {len(data)} values are invalid')

        # Log data statistics
        valid_values = data[valid_data_mask]
        if iFlag_verbose:
            logger.info(f'Data statistics for "{sVariable}":')
            logger.info(f'  - Valid values: {valid_data_count}/{len(data)}')
            logger.info(f'  - Min: {np.min(valid_values):.4f}')
            logger.info(f'  - Max: {np.max(valid_values):.4f}')
            logger.info(f'  - Mean: {np.mean(valid_values):.4f}')
            logger.info(f'  - Std: {np.std(valid_values):.4f}')

        # Create masked connectivity array
        connectivity_masked = np.ma.masked_where(
            self.aConnectivity == -1,
            self.aConnectivity
        )

        # Validate connectivity indices
        valid_connectivity = self.aConnectivity[self.aConnectivity >= 0]
        if len(valid_connectivity) > 0:
            max_vertex_idx = len(self.aVertex_longititude) - 1
            if np.max(valid_connectivity) > max_vertex_idx:
                logger.error(f'Connectivity contains invalid vertex index: max={np.max(valid_connectivity)}, vertices={len(self.aVertex_longititude)}')
                return False

        # Transform to GeoVista unstructured mesh
        if iFlag_verbose:
            logger.info('Creating GeoVista mesh...')
        mesh = gv.Transform.from_unstructured(
            self.aVertex_longititude,
            self.aVertex_latitude,
            connectivity=connectivity_masked,
            crs=crs
        )

        if iFlag_verbose:
            logger.info(f'Created mesh with {mesh.n_cells} cells and {mesh.n_points} points')

        # Attach data to mesh
        scalars = sVariable.capitalize()
        sUnit = sUnit_in if sUnit_in is not None else ""

        # Validate data array length matches mesh cells
        if len(data) != mesh.n_cells:
            logger.error(f'Data array length ({len(data)}) does not match mesh cells ({mesh.n_cells})')
            logger.error(f'This indicates a mismatch between mesh topology and extracted data')
            return False

        mesh.cell_data[scalars] = data
        # Hide/skip cells with NaN by extracting only valid cells

        n_valid = int(np.count_nonzero(valid_data_mask))
        if n_valid == 0:
            logger.warning(f'No valid cells to plot for variable "{scalars}"')
            return False

        valid_cell_indices = np.where(valid_data_mask)[0]
        if iFlag_verbose:
            logger.info(f'Attached data "{scalars}" to mesh cells')

        # Handle animation mode
        if iFlag_create_animation:
            if sFilename_out is None:
                logger.error('Animation mode requires output filename')
                return False

            # Validate animation parameters
            if iAnimation_frames <= 0:
                logger.warning(f'Invalid animation frames {iAnimation_frames}, using default 36')
                iAnimation_frames = 360

            if dAnimation_speed <= 0:
                dAnimation_speed = 360.0 / iAnimation_frames
                if iFlag_verbose:
                    logger.info(f'Calculated animation speed: {dAnimation_speed:.1f}°/frame')

            # Create off-screen plotter for animation
            try:
                plotter = gv.GeoPlotter(off_screen=True)
                if iFlag_verbose:
                    logger.info('Created off-screen plotter for animation')
            except Exception as e:
                logger.error(f'Failed to create off-screen plotter: {e}')
                return False

            # Configure scalar bar
            sargs = {
                "title": f"{scalars} / {sUnit}" if sUnit else scalars,
                "shadow": True,
                "title_font_size": 12,
                "label_font_size": 10,
                "fmt": "%.2f",
                "n_labels": 5,
            }

            # Add mesh to plotter
            mesh_valid = mesh.extract_cells(valid_cell_indices)
            plotter.add_mesh(mesh_valid, scalars=scalars, scalar_bar_args=sargs, cmap=sColormap)

            try:
                # Convert longitude/latitude to radians
                lon_rad = math.radians(dLongitude_focus)
                lat_rad = math.radians(dLatitude_focus)

                # Earth radius (approximately 6371 km, but use normalized units)
                earth_radius = 1.0
                camera_distance = earth_radius * 3.0  # Position camera 3x earth radius away

                # Convert spherical coordinates to Cartesian (x, y, z)
                x_focal = earth_radius * math.cos(lat_rad) * math.cos(lon_rad)
                y_focal = earth_radius * math.cos(lat_rad) * math.sin(lon_rad)
                z_focal = earth_radius * math.sin(lat_rad)

                x_camera = camera_distance * math.cos(lat_rad) * math.cos(lon_rad)
                y_camera = camera_distance * math.cos(lat_rad) * math.sin(lon_rad)
                z_camera = camera_distance * math.sin(lat_rad)

                focal_point = [x_focal, y_focal, z_focal]
                camera_position = [x_camera, y_camera, z_camera]

                plotter.camera.focal_point = focal_point
                plotter.camera.position = camera_position
                plotter.camera.zoom(dZoom_factor)
                if iFlag_verbose:
                    logger.debug(f'Camera configured: focal={focal_point}, position={camera_position}')

            except Exception as e:
                logger.warning(f'Error setting camera position: {e}. Using default view.')

            # Add geographic context
            if iFlag_show_coastlines:
                try:
                    plotter.add_coastlines()
                except Exception as e:
                    logger.warning(f'Could not add coastlines: {e}')

            if iFlag_show_graticule:
                try:
                    plotter.add_graticule(show_labels=True)
                except Exception as e:
                    logger.warning(f'Could not add graticule: {e}')

            # Create animation
            success = _create_rotation_animation(
                plotter, sFilename_out, dLongitude_focus, dLatitude_focus,
                dZoom_factor, iAnimation_frames, dAnimation_speed, sAnimation_format, iFlag_verbose
            )

            plotter.close()
            return success

        else:
            # Single frame visualization
            try:
                if sFilename_out is not None:
                    # Off-screen rendering for saving files
                    plotter = gv.GeoPlotter(off_screen=True)
                    if iFlag_verbose:
                        logger.debug('Created off-screen plotter for file output')
                else:
                    # Interactive rendering
                    plotter = gv.GeoPlotter()
                    if iFlag_verbose:
                        logger.debug('Created interactive plotter')

            except Exception as e:
                logger.error(f'Failed to create GeoVista plotter: {e}')
                return False

            # Configure scalar bar
            sargs = {
                "title": f"{scalars} / {sUnit}" if sUnit else scalars,
                "shadow": True,
                "title_font_size": 12,
                "label_font_size": 10,
                "fmt": "%.2f",
                "n_labels": 5,
            }

            # Add mesh to plotter
            mesh_valid = mesh.extract_cells(valid_cell_indices)
            plotter.add_mesh(mesh_valid, scalars=scalars, scalar_bar_args=sargs, cmap=sColormap)

            # Configure camera position
            try:
                # Convert longitude/latitude to radians
                lon_rad = math.radians(dLongitude_focus)
                lat_rad = math.radians(dLatitude_focus)

                # Earth radius and camera distance
                earth_radius = 1.0
                camera_distance = earth_radius * 3.0

                # Convert spherical coordinates to Cartesian
                x_focal = earth_radius * math.cos(lat_rad) * math.cos(lon_rad)
                y_focal = earth_radius * math.cos(lat_rad) * math.sin(lon_rad)
                z_focal = earth_radius * math.sin(lat_rad)

                x_camera = camera_distance * math.cos(lat_rad) * math.cos(lon_rad)
                y_camera = camera_distance * math.cos(lat_rad) * math.sin(lon_rad)
                z_camera = camera_distance * math.sin(lat_rad)

                focal_point = [x_focal, y_focal, z_focal]
                camera_position = [x_camera, y_camera, z_camera]

                plotter.camera.focal_point = focal_point
                plotter.camera.position = camera_position
                plotter.camera.zoom(dZoom_factor)
                plotter.camera.up = [0, 0, 1]

                if iFlag_verbose:
                    logger.debug(f'Camera configured for focus: ({dLongitude_focus:.2f}°, {dLatitude_focus:.2f}°)')

            except Exception as e:
                logger.warning(f'Error setting camera position: {e}. Using default view.')
                try:
                    plotter.reset_camera()
                except Exception:
                    pass

            # Add geographic context
            if iFlag_show_coastlines:
                try:
                    plotter.add_coastlines()
                    if iFlag_verbose:
                        logger.debug('Added coastlines overlay')
                except Exception as e:
                    logger.warning(f'Could not add coastlines: {e}')

            try:
                plotter.add_axes()
                if iFlag_verbose:
                    logger.debug('Added coordinate axes')
            except Exception as e:
                logger.warning(f'Could not add axes: {e}')

            if iFlag_show_graticule:
                try:
                    plotter.add_graticule(show_labels=True)
                    if iFlag_verbose:
                        logger.debug('Added coordinate graticule with labels')
                except Exception as e:
                    logger.warning(f'Could not add graticule: {e}')

            # Output or display
            if sFilename_out is not None:
                # Save screenshot
                try:
                    plotter.screenshot(sFilename_out)
                    if iFlag_verbose:
                        logger.info(f'✓ Visualization saved to: {sFilename_out}')

                        if os.path.exists(sFilename_out):
                            file_size = os.path.getsize(sFilename_out)
                            logger.info(f'  File size: {file_size / 1024:.1f} KB')
                        else:
                            logger.warning(f'Screenshot command executed but file not found: {sFilename_out}')

                    plotter.close()
                    return True

                except Exception as e:
                    logger.error(f'Failed to save screenshot: {e}')
                    logger.error(f'Traceback: {traceback.format_exc()}')
                    plotter.close()
                    return False
            else:
                # Interactive display
                try:
                    if iFlag_verbose:
                        logger.info('Opening interactive visualization window...')
                    plotter.show()
                    return True
                except Exception as e:
                    logger.error(f'Failed to display interactive visualization: {e}')
                    logger.error(f'Traceback: {traceback.format_exc()}')
                    plotter.close()
                    return False

    except ImportError as e:
        logger.error(f'Missing required GeoVista dependencies: {e}')
        return False

    except Exception as e:
        logger.error(f'Unexpected error during target mesh visualization: {e}')
        logger.error(f'Error type: {type(e).__name__}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return False


def _create_rotation_animation(plotter, sFilename_out, dLongitude_start, dLatitude_focus,
                               dZoom_factor, iAnimation_frames, dAnimation_speed, sAnimation_format, iFlag_verbose):
    """
    Create a rotating animation of the 3D globe visualization.

    Generates multiple frames by rotating the camera around the globe at a fixed latitude,
    then combines them into a video file using imageio.

    Args:
        plotter: GeoVista plotter instance with mesh already added
        sFilename_out (str): Output animation file path (e.g., 'animation.mp4')
        dLongitude_start (float): Starting longitude for rotation in degrees
        dLatitude_focus (float): Fixed latitude for camera focus in degrees
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
            logger.info(f'  - Fixed latitude: {dLatitude_focus:.1f}°')
            logger.info(f'  - Rotation speed: {dAnimation_speed:.1f}°/frame')
            logger.info(f'  - Output format: {sAnimation_format}')


        earth_radius = 1.0
        camera_distance = earth_radius * 3.0
        # Generate frames
        amplitude_deg = 20.0
        cycles = 1.0
        phase = 0.0
        plotter.open_movie(sFilename_out, framerate=30)
        for iFrame in range(iAnimation_frames):
            # Calculate current longitude
            dLongitude_current = dLongitude_start + (iFrame * dAnimation_speed)
            dLongitude_current = dLongitude_current % 360.0  # Keep within [0, 360)
            if dLongitude_current > 180.0:
                dLongitude_current -= 360.0  # Convert to [-180, 180]

            # Sine-pattern latitude: oscillate latitude as a sine curve while rotating longitude.
            # amplitude (deg) controls north/south swing; cycles controls how many sine cycles
            # occur over the full animation. phase shifts the starting point if needed.
            # Protect against division by zero
            frames_div = float(iAnimation_frames) if iAnimation_frames > 0 else 1.0
            theta = 2.0 * math.pi * (float(iFrame) / frames_div) * cycles + phase
            dLatitude_current = float(dLatitude_focus) + amplitude_deg * math.sin(theta)
            # Clip to avoid exactly hitting poles (which can cause view issues)
            dLatitude_current = max(-89.9, min(89.9, dLatitude_current))

            # Convert to radians
            lon_rad = math.radians(dLongitude_current)
            lat_rad = math.radians(dLatitude_current)

            # Calculate camera position
            x_focal = earth_radius * math.cos(lat_rad) * math.cos(lon_rad)
            y_focal = earth_radius * math.cos(lat_rad) * math.sin(lon_rad)
            z_focal = earth_radius * math.sin(lat_rad)

            x_camera = camera_distance * math.cos(lat_rad) * math.cos(lon_rad)
            y_camera = camera_distance * math.cos(lat_rad) * math.sin(lon_rad)
            z_camera = camera_distance * math.sin(lat_rad)

            focal_point = [x_focal, y_focal, z_focal]
            camera_position = [x_camera, y_camera, z_camera]

            # Update camera
            plotter.camera.focal_point = focal_point
            plotter.camera.position = camera_position
            #plotter.camera.zoom(dZoom_factor)
            plotter.camera.up = [0, 0, 1]
            plotter.add_axes()  # Re-add axes to ensure visibility
            plotter.render()  # Render the scene

            # Render frame
            try:
                #use write_frame to get the image as a numpy array
                plotter.write_frame()

            except Exception as e:
                logger.error(f'Failed to render frame {iFrame + 1}: {e}')
                return False

        # Determine output format and settings
        if sAnimation_format.lower() == 'gif':
            pass

        elif sAnimation_format.lower() in ['mp4', 'avi']:
            # Video format
            output_file = sFilename_out
            if sAnimation_format.lower() == 'mp4' and not output_file.lower().endswith('.mp4'):
                output_file = os.path.splitext(sFilename_out)[0] + '.mp4'
            elif sAnimation_format.lower() == 'avi' and not output_file.lower().endswith('.avi'):
                output_file = os.path.splitext(sFilename_out)[0] + '.avi'

            try:
                if iFlag_verbose:
                    logger.info(f'✓ Video animation saved to: {output_file}')
            except Exception as e:
                logger.error(f'Failed to save video animation: {e}')
                logger.error('Ensure ffmpeg is installed for MP4 support: pip install imageio[ffmpeg]')
                return False

        else:
            logger.error(f'Unsupported animation format: {sAnimation_format}')
            logger.error('Supported formats: mp4, gif, avi')
            return False

        if os.path.exists(sFilename_out):
            file_size = os.path.getsize(sFilename_out)
            if iFlag_verbose:
                logger.info(f'✓ Animation created successfully: {sFilename_out}')
                logger.info(f'  File size: {file_size / (1024*1024):.2f} MB')
                logger.info(f'  Frames: {iAnimation_frames}')
                logger.info(f'  Format: {sAnimation_format.upper()}')
            return True
        else:
            logger.error('Animation file was not created')
            return False

    except Exception as e:
        logger.error(f'Unexpected error during animation creation: {e}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return False
