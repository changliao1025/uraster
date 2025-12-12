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
gdal.UseExceptions()
from pyearth.visual.geovista.utility import (
    VisualizationConfig,
    AnimationConfig,
    setup_geovista_plotter,
    configure_camera,
    add_geographic_context,
)
from pyearth.visual.geovista.map_single_frame import map_single_frame
from pyearth.visual.geovista.animate_rotating_frames import animate_rotating_frames
from uraster.classes.sraster import sraster
from uraster import utility
from uraster.utility import setup_logger
logger = setup_logger(__name__.split('.')[-1])
CRS = "EPSG:4326"

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
                    iFlag_verbose_in: bool = False) -> bool:
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
        iFlag_verbose_in: If True, print detailed progress messages.

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
        verbose=iFlag_verbose_in
    )

    try:
        import geovista as gv

        if config.verbose:
            logger.info(f'Visualizing {len(self.aFilename_source_raster)} raster file(s)...')

        # Setup plotter


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



        # Handle output
        map_single_frame(pMesh_raster, sScalar, aValid_cell_indices ,sUnit, config, sFilename_out)

    except ImportError as e:
        logger.error('GeoVista library not available. Install with: pip install geovista')
        logger.error(f'Import error: {e}')
        return False

    except Exception as e:
        logger.error(f'Unexpected error during raster visualization: {e}')
        logger.error(f'Traceback: {traceback.format_exc()}')
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
                          iFlag_wireframe_only: bool = False,
                          dEdge_width: float = 1.0,
                          sEdge_color: str = 'black',
                          iFlag_verbose_in: bool = False) -> bool:
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
        iFlag_wireframe_only: If True, show only mesh edges without cell filling.
            Default is False (shows filled cells).
        dEdge_width: Line width for mesh edges when wireframe mode is enabled.
            Default is 1.0.
        sEdge_color: Color for mesh edges in wireframe mode. Default is 'black'.
            Examples: 'white', 'red', 'blue', 'gray', or RGB tuples like (1.0, 0.0, 0.0).
        iFlag_verbose_in: If True, print detailed progress messages. Default is False.

    Returns:
        True if visualization successful, False otherwise

    Note:
        - Requires 'geovista' package: pip install geovista
        - Interactive mode requires display environment
        - Mesh topology must be built before visualization (call rebuild_mesh_topology first)
        - Wireframe mode is useful for examining mesh structure and topology
    """
    # Validate inputs using new utility functions
    if not _validate_mesh_data(self):
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
        verbose=iFlag_verbose_in
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
        pMesh_source = gv.Transform.from_unstructured(
            self.aVertex_longititude,
            self.aVertex_latitude,
            connectivity=connectivity_masked,
            crs=CRS
        )

        # Validate cell data array length matches mesh cells
        if len(self.aCellID) != pMesh_source.n_cells:
            logger.error(f'Cell ID array length ({len(self.aCellID)}) does not match '
                        f'mesh cells ({pMesh_source.n_cells})')
            return False

        # Prepare mesh metadata
        sScalar = 'Cell ID'
        pMesh_source.cell_data[sScalar] = self.aCellID

        if config.verbose:
            logger.info(f'Created GeoVista mesh with {pMesh_source.n_cells} cells and {pMesh_source.n_points} points')

        # Get valid cell indices (non-NaN values)
        aValid_data_mask = np.isfinite(self.aCellID)
        iFlag_valid = int(np.count_nonzero(aValid_data_mask))

        if iFlag_valid == 0:
            logger.warning(f'No valid cells to plot for variable "{sScalar}"')
            return None

        aValid_cell_indices = np.where(aValid_data_mask)[0]
        sUnit=''
        map_single_frame(pMesh_source, sScalar, aValid_cell_indices , config,
                         sUnit=sFilename_out, sFilename_out = sFilename_out)

        # Output or display
        return

    except ImportError as e:
        logger.error('GeoVista library not available. Install with: pip install geovista')
        logger.error(f'Import error: {e}')
        return False

    except Exception as e:
        logger.error(f'Unexpected error during mesh visualization: {e}')
        logger.error(f'Error type: {type(e).__name__}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return False

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
                         iFlag_verbose_in: Optional[bool] = False) -> bool:
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
        iFlag_verbose_in (bool, optional): If True, print detailed progress messages.
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
        verbose=iFlag_verbose_in
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

        aData, nFeature = pData_result


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
            animate_rotating_frames(pMesh, sScalars, aValid_cell_indices, sUnit, config, sFilename_out)
        else:

            map_single_frame(pMesh, sScalars, aValid_cell_indices, sUnit, config, sFilename_out)

    except ImportError as e:
        logger.error('GeoVista library not available. Install with: pip install geovista')
        logger.error(f'Import error: {e}')
        return False

    except Exception as e:
        logger.error(f'Unexpected error during target mesh visualization: {e}')
        logger.error(f'Error type: {type(e).__name__}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return False

def _validate_mesh_data(pUraster_instance) -> bool:
    """
    Validate that mesh data is available and properly formatted.

    Args:
        pUraster_instance: The uraster instance to validate

    Returns:
        bool: True if mesh data is valid, False otherwise
    """
    if pUraster_instance.aVertex_longititude is None or pUraster_instance.aVertex_latitude is None:
        logger.error('Mesh vertices not available. Build mesh topology first.')
        return False

    if pUraster_instance.aConnectivity is None:
        logger.error('Mesh connectivity not available. Build mesh topology first.')
        return False

    if len(pUraster_instance.aVertex_longititude) == 0 or len(pUraster_instance.aVertex_latitude) == 0:
        logger.error('Mesh vertices are empty.')
        return False

    if pUraster_instance.aConnectivity.size == 0:
        logger.error('Mesh connectivity is empty.')
        return False

    return True

def _extract_target_mesh_data(sFilename: str, sVariable: str, iFlag_verbose_in: bool = False) -> Optional[Tuple[np.ndarray, int]]:
    """
    Extract variable data from target mesh file.

    Args:
        sFilename: Path to target mesh file
        sVariable: Variable name to extract
        iFlag_verbose_in: Enable verbose logging

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

        if iFlag_verbose_in:
            logger.info(f'Target mesh contains {nFeatures} features with {iFieldCount} fields')

        # Check if variable field exists
        aField_names = [pLayerDefn.GetFieldDefn(i).GetName() for i in range(iFieldCount)]
        if sVariable not in aField_names:
            logger.error(f'Variable "{sVariable}" not found in target mesh')
            logger.error(f'Available fields: {", ".join(aField_names)}')
            pDataset = None
            return None

        if iFlag_verbose_in:
            logger.info(f'Extracting variable: {sVariable}')

        # Extract variable data from features, handling multipolygons correctly
        aData_list = []
        pLayer.ResetReading()
        pFeature = pLayer.GetNextFeature()
        iFeature_count = 0

        iCount_multipolygons = 0

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
                    iCount_multipolygons += 1
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
                            else:
                                pPolygon_part.FlattenTo2D()
                                sWkt = pPolygon_part.ExportToWkt()
                                logger.warning(f'Invalid polygon part {iPart} in multipolygon feature {iFeature_count}: {sWkt} ')
                                aData_list.append(np.nan)

                        if iValid_parts == 0:
                            logger.warning(f'No valid parts found in multipolygon feature {iFeature_count}')
                            #aData_list.append(np.nan)

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

def _create_target_mesh(pUraster_instance, aData: np.ndarray, sVariable: str, iFlag_verbose_in: bool = False) -> Optional[Tuple[Any, str, np.ndarray]]:
    """
    Create GeoVista mesh from uraster instance and attach data.

    Args:
        pUraster_instance: The uraster instance
        aData: Data array to attach to mesh
        sVariable: Variable name for the data
        iFlag_verbose_in: Enable verbose logging

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
        if iFlag_verbose_in:
            logger.info('Creating GeoVista mesh...')

        pMesh = gv.Transform.from_unstructured(
            pUraster_instance.aVertex_longititude,
            pUraster_instance.aVertex_latitude,
            connectivity=aConnectivity_masked,
            crs=CRS
        )

        if iFlag_verbose_in:
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

        if iFlag_verbose_in:
            logger.info(f'Attached data "{sScalars}" to mesh cells')
            logger.info(f'Valid cells for visualization: {iN_valid}/{len(aData)}')

        return pMesh, sScalars, aValid_cell_indices

    except Exception as e:
        logger.error(f'Error creating target mesh: {e}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return None




