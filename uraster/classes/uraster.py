import os
import logging
import traceback
import numpy as np
from osgeo import gdal, ogr, osr
gdal.UseExceptions()
from uraster.operation import extract, intersect
from uraster.classes import _visual
from uraster.classes.sraster import sraster
# Set up logging for crash detection
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)
crs = "EPSG:4326"
pDriver_geojson = ogr.GetDriverByName('GeoJSON')
pDriver_shp = ogr.GetDriverByName('ESRI Shapefile')

class uraster:
    """
    Unstructured raster processing class for zonal statistics on mesh geometries.

    Handles complex scenarios including:
    - International Date Line (IDL) crossing polygons
    - Parallel processing for large datasets
    - Multiple raster formats and coordinate systems
    - Comprehensive error handling and crash detection
    """

    def __init__(self, aConfig=None):
        """
        Initialize uraster instance.

        Args:
            aConfig (dict, optional): Configuration dictionary with keys:
                - iFlag_remap_method (int): Remap method (1=nearest, 2=nearest, 3=weighted average)
                - sFilename_source_mesh (str): Source mesh file path
                - sFilename_target_mesh (str): Target mesh file path
                - aFilename_source_raster (list): List of source raster file paths
        """
        # Default configuration
        if aConfig is None:
            aConfig = {}

        # Processing flags and resolutions
        self.iFlag_global = None
        self.iFlag_remap_method = aConfig.get(
            'iFlag_remap_method', 1)  # Default to nearest neighbor
        self.dResolution_raster = None
        self.dResolution_uraster = None

        # File paths
        self.sFilename_source_mesh = aConfig.get('sFilename_source_mesh', None)
        self.sField_unique_id = aConfig.get('sField_unique_id', None)
        self.sFilename_target_mesh = aConfig.get('sFilename_target_mesh', None)
        self.aFilename_source_raster = aConfig.get(
            'aFilename_source_raster', [])

        # Cell counts
        self.nCell = -1
        self.nCell_source = -1
        self.nCell_target = -1
        self.nVertex_max = 0  # Will be calculated dynamically

        # Mesh topology data
        self.aVertex_longititude = None
        self.aVertex_latitude = None
        self.aCenter_longititude = None
        self.aCenter_latitude = None
        self.aConnectivity = None
        self.aCellID = None

        # Mesh area statistics
        self.dArea_min = None
        self.dArea_max = None
        self.dArea_mean = None

        # Resolution comparison threshold (ratio of mesh to raster resolution)
        # If mesh cells are within this factor of raster resolution, use weighted averaging
        # mesh resolution < 3x raster resolution triggers weighted avg
        self.dResolution_ratio_threshold = 3.0

        # Validate configuration
        if self.iFlag_remap_method not in [1, 2, 3]:
            logger.warning(
                f"Invalid remap method {self.iFlag_remap_method}, defaulting to 1 (nearest neighbor)")
            self.iFlag_remap_method = 1

    def setup(self, iFlag_verbose=False):
        """
        Initialize and validate the uraster configuration.
        Checks raster files and mesh file for existence and validity.

        Args:
            iFlag_verbose (bool, optional): If True, print detailed progress messages.
                If False, only print error messages. Default is False.

        Returns:
            bool: True if setup successful, False otherwise
        """
        raster_check = self.check_raster_files(iFlag_verbose=iFlag_verbose)
        mesh_check = self.check_mesh_file(iFlag_verbose=iFlag_verbose)

        return raster_check is not None and mesh_check is not None

    def check_raster_files(self, aFilename_source_raster_in=None, iFlag_verbose=False):
        """
        Validate and prepare input raster files, converting to WGS84 if needed.

        Performs comprehensive validation of raster files including:
        - File existence and readability
        - Valid GDAL raster format
        - Coordinate system compatibility
        - Data integrity checks

        Args:
            aFilename_source_raster_in (list, optional): List of raster file paths.
                If None, uses self.aFilename_source_raster
            iFlag_verbose (bool, optional): If True, print detailed progress messages.
                If False, only print error messages. Default is False.

        Returns:
            list: List of WGS84 raster file paths, or None if validation fails

        Note:
            - Non-WGS84 rasters are automatically converted and cached
            - All rasters must be valid and readable for processing to continue
        """
        # Determine input raster list
        if aFilename_source_raster_in is None:
            aFilename_source_raster = self.aFilename_source_raster
        else:
            aFilename_source_raster = aFilename_source_raster_in

        # Validate input list
        if not aFilename_source_raster:
            logger.error('No raster files provided for validation')
            return None

        if not isinstance(aFilename_source_raster, (list, tuple)):
            logger.error(
                f'Raster files must be provided as a list, got {type(aFilename_source_raster).__name__}')
            return None

        if iFlag_verbose:
            logger.info(
                f'Validating {len(aFilename_source_raster)} raster file(s)...')

        # Phase 1: Check file existence and readability
        for idx, sFilename_raster_in in enumerate(aFilename_source_raster, 1):
            if not isinstance(sFilename_raster_in, str):
                logger.error(
                    f'Raster file path must be a string, got {type(sFilename_raster_in).__name__} at index {idx}')
                return None

            if not sFilename_raster_in.strip():
                logger.error(f'Empty raster file path at index {idx}')
                return None

            if not os.path.exists(sFilename_raster_in):
                logger.error(
                    f'Raster file does not exist: {sFilename_raster_in}')
                return None

            if not os.path.isfile(sFilename_raster_in):
                logger.error(f'Path is not a file: {sFilename_raster_in}')
                return None

            # Check file permissions
            if not os.access(sFilename_raster_in, os.R_OK):
                logger.error(
                    f'Raster file is not readable: {sFilename_raster_in}')
                return None

            # Quick GDAL format validation
            try:
                pDataset_test = gdal.Open(
                    sFilename_raster_in, gdal.GA_ReadOnly)
                if pDataset_test is None:
                    logger.error(
                        f'GDAL cannot open raster file: {sFilename_raster_in}')
                    return None
                pDataset_test = None  # Close dataset
            except Exception as e:
                logger.error(
                    f'Error opening raster with GDAL: {sFilename_raster_in}: {e}')
                return None

        if iFlag_verbose:
            logger.info('All raster files exist and are readable')

        # Phase 2: Process and convert rasters to WGS84
        aFilename_source_raster_out = []

        # Create WGS84 spatial reference for comparison
        pSpatialRef_wgs84 = None
        try:
            pSpatialRef_wgs84 = osr.SpatialReference()
            pSpatialRef_wgs84.ImportFromEPSG(4326)
            wkt_wgs84 = pSpatialRef_wgs84.ExportToWkt()
        except Exception as e:
            logger.error(f'Failed to create WGS84 spatial reference: {e}')
            return None
        finally:
            # Clean up spatial reference object
            if pSpatialRef_wgs84 is not None:
                pSpatialRef_wgs84 = None

        # Process each raster file
        for idx, sFilename_raster_in in enumerate(aFilename_source_raster, 1):
            if iFlag_verbose:
                logger.info(
                    f'Processing raster {idx}/{len(aFilename_source_raster)}: {os.path.basename(sFilename_raster_in)}')

            try:
                # Create sraster instance and read metadata
                pRaster = sraster(sFilename_raster_in)
                pRaster.read_metadata()

                # Validate critical metadata
                if pRaster.pSpatialRef_wkt is None:
                    logger.error(
                        f'Raster has no spatial reference: {sFilename_raster_in}')
                    return None

                if pRaster.nrow is None or pRaster.ncolumn is None:
                    logger.error(
                        f'Invalid raster dimensions: {sFilename_raster_in}')
                    return None

                if pRaster.nrow <= 0 or pRaster.ncolumn <= 0:
                    logger.error(
                        f'Raster has invalid dimensions ({pRaster.nrow}x{pRaster.ncolumn}): {sFilename_raster_in}')
                    return None

                # Check if coordinate system matches WGS84
                if pRaster.pSpatialRef_wkt == wkt_wgs84:
                    if iFlag_verbose:
                        logger.info(f'  ✓ Already in WGS84 (EPSG:4326)')
                    aFilename_source_raster_out.append(sFilename_raster_in)
                else:
                    # Convert to WGS84
                    if iFlag_verbose:
                        logger.info(
                            f'  → Converting to WGS84 from {pRaster.pSpatialRef.GetName() if pRaster.pSpatialRef else "unknown CRS"}')
                    try:
                        pRaster_wgs84 = pRaster.convert_to_wgs84()

                        if pRaster_wgs84 is None or not hasattr(pRaster_wgs84, 'sFilename'):
                            logger.error(
                                f'Conversion to WGS84 failed: {sFilename_raster_in}')
                            return None

                        if not os.path.exists(pRaster_wgs84.sFilename):
                            logger.error(
                                f'Converted WGS84 file not found: {pRaster_wgs84.sFilename}')
                            return None

                        if iFlag_verbose:
                            logger.info(
                                f'  ✓ Converted to: {pRaster_wgs84.sFilename}')
                        aFilename_source_raster_out.append(
                            pRaster_wgs84.sFilename)

                    except Exception as e:
                        logger.error(
                            f'Error during WGS84 conversion: {sFilename_raster_in}: {e}')
                        logger.error(f'Traceback: {traceback.format_exc()}')
                        return None

                # Log raster summary
                if iFlag_verbose:
                    logger.debug(
                        f'  - Dimensions: {pRaster.nrow} x {pRaster.ncolumn} pixels')
                    logger.debug(f'  - Data type: {pRaster.eType}')
                    if hasattr(pRaster, 'dNoData'):
                        logger.debug(f'  - NoData value: {pRaster.dNoData}')

                pRaster.pSpatialRef = None  # Clean up spatial reference

            except AttributeError as e:
                logger.error(
                    f'Missing expected attribute in sraster: {sFilename_raster_in}: {e}')
                logger.error(
                    f'Ensure sraster class has all required methods and attributes')
                return None

            except Exception as e:
                logger.error(
                    f'Unexpected error processing raster {sFilename_raster_in}: {e}')
                logger.error(f'Error type: {type(e).__name__}')
                logger.error(f'Traceback: {traceback.format_exc()}')
                return None

        # Final validation
        if len(aFilename_source_raster_out) != len(aFilename_source_raster):
            logger.error(
                f'Output count mismatch: expected {len(aFilename_source_raster)}, got {len(aFilename_source_raster_out)}')
            return None

        if iFlag_verbose:
            logger.info(
                f'Successfully validated and prepared {len(aFilename_source_raster_out)} raster file(s)')

        return aFilename_source_raster_out

    def check_mesh_file(self, iFlag_verbose=False):
        """
        Check if the source mesh file exists and build its topology.

        Args:
            iFlag_verbose (bool, optional): If True, print detailed progress messages.
                If False, only print error messages. Default is False.

        Returns:
            tuple or None: (vertices_lon, vertices_lat, connectivity) if successful, None otherwise
        """
        if not self.sFilename_source_mesh:
            logger.error("No source mesh filename provided")
            return None

        if not os.path.exists(self.sFilename_source_mesh):
            logger.error(
                f"Source mesh file does not exist: {self.sFilename_source_mesh}")
            return None

        return self.rebuild_mesh_topology(iFlag_verbose=iFlag_verbose)

    def _get_geometry_type_name(self, geometry_type):
        """
        Convert OGR geometry type integer to readable string name.

        Handles both standard and 3D/Z-flagged geometry types.

        Args:
            geometry_type (int): OGR geometry type constant

        Returns:
            str: Human-readable geometry type name (e.g., "wkbPolygon")
        """
        geometry_types = {
            ogr.wkbUnknown: "wkbUnknown",
            ogr.wkbPoint: "wkbPoint",
            ogr.wkbLineString: "wkbLineString",
            ogr.wkbPolygon: "wkbPolygon",
            ogr.wkbMultiPoint: "wkbMultiPoint",
            ogr.wkbMultiLineString: "wkbMultiLineString",
            ogr.wkbMultiPolygon: "wkbMultiPolygon",
            ogr.wkbGeometryCollection: "wkbGeometryCollection"
        }

        # Direct match
        if geometry_type in geometry_types:
            return geometry_types[geometry_type]

        # Check base type (removes 3D/Z flags)
        base_type = geometry_type & 0xFF
        for const_val, name in geometry_types.items():
            if (const_val & 0xFF) == base_type:
                return f"{name} (with flags)"

        return f"Unknown geometry type: {geometry_type}"

    def rebuild_mesh_topology(self, iFlag_verbose=False):
        """
        Rebuild mesh topology from source mesh file by extracting vertices,
        connectivity, and centroids for unstructured mesh processing.

        This method uses the enhanced standalone rebuild_mesh_topology function
        and updates all instance attributes with the comprehensive mesh information.

        Args:
            iFlag_verbose (bool, optional): If True, print detailed progress messages.
                If False, only print error messages. Default is False.

        Returns:
            tuple: (vertices_longitude, vertices_latitude, connectivity) or None on failure
        """
        # Use the enhanced standalone function from _visual module
        mesh_info = _visual.rebuild_mesh_topology(
            self.sFilename_source_mesh,
            iFlag_verbose=iFlag_verbose,
            sField_unique_id=self.sField_unique_id
        )

        if mesh_info is None:
            return None

        # Update all instance attributes with comprehensive mesh information
        self.aVertex_longititude = mesh_info['vertices_longitude']
        self.aVertex_latitude = mesh_info['vertices_latitude']
        self.aConnectivity = mesh_info['connectivity']
        self.aCenter_longititude = mesh_info['cell_centroids_longitude']
        self.aCenter_latitude = mesh_info['cell_centroids_latitude']
        self.aCellID = mesh_info['cell_ids']
        self.dArea_min = mesh_info['area_min']
        self.dArea_max = mesh_info['area_max']
        self.dArea_mean = mesh_info['area_mean']
        self.nVertex_max = mesh_info['max_vertices_per_cell']
        self.nCell_source = mesh_info['num_cells']  # Update cell count

        # Return the traditional tuple format for backward compatibility
        return (mesh_info['vertices_longitude'],
                mesh_info['vertices_latitude'],
                mesh_info['connectivity'])

    def report_inputs(self, iFlag_show_gpu_info=False):
        """
        Print comprehensive input information including raster and mesh details.

        Args:
            iFlag_show_gpu_info (bool): If True, also print GPU/GeoVista information
        """
        self.print_raster_info()
        self.print_mesh_info()

        if iFlag_show_gpu_info:
            try:
                import geovista.report as gvreport
                print("\n" + "="*60)
                print("GPU/GeoVista Information:")
                print("="*60)
                print(gvreport.Report())
            except ImportError:
                logger.warning("GeoVista not available for GPU info reporting")

    def report_outputs(self, sFilename_output=None):
        """
        Report output statistics.

        Args:
            sFilename_output (str, optional): Output file to report on
        """
        if sFilename_output and os.path.exists(sFilename_output):
            logger.info(f"Output file created: {sFilename_output}")
            logger.info(
                f"Output file size: {os.path.getsize(sFilename_output) / (1024*1024):.2f} MB")
        else:
            logger.warning("No output file information available")

    def print_raster_info(self):
        """
        Print detailed information about all input raster files.
        """
        print("\n" + "="*60)
        print(
            f"Input Raster Information ({len(self.aFilename_source_raster)} file(s)):")
        print("="*60)

        for idx, sFilename in enumerate(self.aFilename_source_raster, 1):
            print(f"\n[{idx}] {sFilename}")
            try:
                pRaster = sraster(sFilename)
                pRaster.read_metadata()
                pRaster.print_info()
            except Exception as e:
                logger.error(f"Error reading raster info: {e}")

    def print_mesh_info(self):
        """
        Print detailed mesh topology information.
        """
        if self.aCenter_longititude is None or len(self.aCenter_longititude) == 0:
            logger.warning("Mesh topology not yet built")
            return

        print("\n" + "="*60)
        print("Mesh Topology Information:")
        print("="*60)
        print(f"Number of mesh cells: {len(self.aCenter_longititude)}")
        print(
            f"Cell longitude range: {self.aCenter_longititude.min():.3f} to {self.aCenter_longititude.max():.3f}")
        print(
            f"Cell latitude range: {self.aCenter_latitude.min():.3f} to {self.aCenter_latitude.max():.3f}")
        print(f"Maximum vertices per cell: {self.nVertex_max}")

        if self.aVertex_longititude is not None:
            print(f"Total unique vertices: {len(self.aVertex_longititude)}")
        if self.aConnectivity is not None:
            print(f"Connectivity matrix shape: {self.aConnectivity.shape}")

        # Display area statistics if available
        if self.dArea_min is not None and self.dArea_max is not None:
            print(f"\nCell Area Statistics:")
            print(f"  Min area: {self.dArea_min:.6f}")
            print(f"  Max area: {self.dArea_max:.6f}")
            print(f"  Mean area: {self.dArea_mean:.6f}")
            print(
                f"  Area range ratio: {self.dArea_max/self.dArea_min:.2f}x" if self.dArea_min > 0 else "  Area range ratio: N/A")

        print("="*60)

    def run_remap(self, sFilename_target_mesh_out=None,
                  sFilename_source_mesh_in=None,
                  aFilename_source_raster_in=None,
                  iFlag_stat_in=True,
                  iFlag_weighted_average_in=False,
                  iFlag_remap_method_in=1,
                  iFlag_save_clipped_raster_in=0,
                  sFolder_raster_out_in=None,
                  iFlag_discrete_in=False,
                  iFlag_verbose=False):
        """
        Perform zonal statistics by clipping raster data to mesh polygons.

        This method delegates to the extract module for implementation.
        """
        if aFilename_source_raster_in is None:
            aFilename_source_raster = self.aFilename_source_raster
        else:
            aFilename_source_raster = aFilename_source_raster_in

        if sFilename_source_mesh_in is None:
            sFilename_source_mesh = self.sFilename_source_mesh
        else:
            sFilename_source_mesh = sFilename_source_mesh_in

        if sFilename_target_mesh_out is None:
            sFilename_target_mesh = self.sFilename_target_mesh
        else:
            sFilename_target_mesh = sFilename_target_mesh_out
            self.sFilename_target_mesh = sFilename_target_mesh_out


        #check stat and discrete compatibility
        if iFlag_discrete_in:
            #for discrete, only remap method 1 (nearest neighbor) is allowed
            #can we apply statistics for discrete?
            iFlag_stat_in = False
            if iFlag_remap_method_in != 1:
                logger.error("For discrete remap, only remap method 1 (nearest neighbor) is allowed.")
                return None
        else:
            #for continuous, all remap methods are allowed
            pass

        # the model should suport weighted average and discrete remap
        if iFlag_weighted_average_in:
            # call the polygon calculation with weighted average
            sFilename_raster = aFilename_source_raster[0]
            pRaster = sraster(sFilename_raster)
            pRaster.read_metadata()
            sFilename_raster_mesh = pRaster.create_raster_mesh()
            return intersect.run_remap(
                sFilename_target_mesh,
                sFilename_source_mesh,
                sFilename_raster,
                sFilename_raster_mesh,
                self.dArea_min,
                iFlag_save_clipped_raster_in=iFlag_save_clipped_raster_in,
                sFolder_raster_out_in=sFolder_raster_out_in,
                iFlag_discrete_in=iFlag_discrete_in,
                iFlag_verbose=iFlag_verbose)

        else:
            return extract.run_remap(
                sFilename_target_mesh,
                sFilename_source_mesh,
                aFilename_source_raster,
                self.dArea_min,
                iFlag_remap_method_in=iFlag_remap_method_in,
                iFlag_discrete_in=iFlag_discrete_in,
                iFlag_stat_in=iFlag_stat_in,
                iFlag_save_clipped_raster_in=iFlag_save_clipped_raster_in,
                sFolder_raster_out_in=sFolder_raster_out_in,
                iFlag_verbose_in=iFlag_verbose)

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
        return _visual.visualize_source_mesh(
            self, sFilename_out, dLongitude_focus_in=dLongitude_focus_in, dLatitude_focus_in=dLatitude_focus_in,
            dZoom_factor=dZoom_factor,
            iFlag_show_coastlines=iFlag_show_coastlines,
            iFlag_show_graticule=iFlag_show_graticule,
            iFlag_verbose=iFlag_verbose
        )

    def visualize_raster(self, sFilename_out=None, iFlag_verbose=False):
        """
        Visualize source raster data using GeoVista.

        Note:
            Not yet implemented. Placeholder for future raster visualization.
        """
        return _visual.visualize_raster(self, sFilename_out=sFilename_out, iFlag_verbose=iFlag_verbose)

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
                              iAnimation_frames=360,
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
                Default is 10.0. Calculated as 360 / iAnimation_frames if not specified.
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
        return _visual.visualize_target_mesh(
            self, sVariable_in, sUnit_in, sFilename_out,
            dLongitude_focus_in=dLongitude_focus_in,
            dLatitude_focus_in=dLatitude_focus_in,
            dZoom_factor=dZoom_factor,
            iFlag_show_coastlines=iFlag_show_coastlines,
            iFlag_show_graticule=iFlag_show_graticule,
            sColormap=sColormap,
            iFlag_create_animation=iFlag_create_animation,
            iAnimation_frames=iAnimation_frames,
            dAnimation_speed=dAnimation_speed,
            sAnimation_format=sAnimation_format,
            iFlag_verbose=iFlag_verbose
        )

    def _create_rotation_animation(self, plotter, sFilename_out, dLongitude_start, dLatitude_focus,
                                   iAnimation_frames, dAnimation_speed, sAnimation_format, iFlag_verbose=False):
        """
        Create a rotating animation of the 3D globe visualization.

        This method delegates to the _visual module for implementation.
        """
        return _visual._create_rotation_animation(
            self, plotter, sFilename_out, dLongitude_start, dLatitude_focus,
            iAnimation_frames, dAnimation_speed, sAnimation_format, iFlag_verbose
        )

    def cleanup(self):
        """
        Cleanup method to release spatial reference objects and other resources.
        """
        try:
            if hasattr(self, 'pSpatialRef') and self.pSpatialRef is not None:
                self.pSpatialRef = None
                logger.debug(
                    'Spatial reference object cleaned up successfully')
        except Exception as e:
            logger.warning(f'Error during cleanup of spatial reference: {e}')
