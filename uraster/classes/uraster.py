import os
import logging
import traceback
import numpy as np
from osgeo import gdal, ogr, osr

gdal.UseExceptions()
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename
from uraster.operation import extract, intersect
from uraster.classes import _visual
from uraster.classes.sraster import sraster
from uraster import utility
from uraster.utility import setup_logger

logger = setup_logger(__name__.split(".")[-1])
crs = "EPSG:4326"
pDriver_geojson = ogr.GetDriverByName("GeoJSON")
pDriver_shp = ogr.GetDriverByName("ESRI Shapefile")


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
            "iFlag_remap_method", 1
        )  # Default to nearest neighbor
        self.dResolution_raster = None
        self.dResolution_uraster = None

        # File paths
        self.sFilename_source_mesh = aConfig.get("sFilename_source_mesh", None)
        self.sField_unique_id = aConfig.get(
            "sField_unique_id", "cellid"
        )  # Default to 'cellid'
        self.iField_unique_type = None  # Will be set during setup_mesh_cellid
        self.sFilename_target_mesh = aConfig.get("sFilename_target_mesh", None)
        self.aFilename_source_raster = aConfig.get("aFilename_source_raster", [])

        self.iFlag_polar = aConfig.get("iFlag_polar", 0)  # Default to 0 (non-polar)
        self.iFlag_global = aConfig.get("iFlag_global", 1)  # Default to 1 (global)
        self.iFlag_discrete = aConfig.get(
            "iFlag_discrete", 0
        )  # Default to 0 (continuous)

        # Cell counts
        self.nCell = -1
        self.nCell_source = -1
        self.nCell_target = -1
        self.nPolygon = -1
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
                f"Invalid remap method {self.iFlag_remap_method}, defaulting to 1 (nearest neighbor)"
            )
            self.iFlag_remap_method = 1

        # Memory management: raster caching infrastructure
        self._raster_cache = {}  # filename -> sraster instance (metadata only)
        self._cache_size_threshold = (
            100 * 1024 * 1024
        )  # 100MB threshold for caching decisions
        self._cache_enabled = True  # Allow disabling cache if needed

    def _get_sraster(
        self,
        sFilename_or_aFilename,
        load_data=False,
        use_cache=True,
        iFlag_verbose_in=False,
    ):
        """
        Get sraster instance(s) with intelligent caching for single files or tiled datasets.

        Args:
            sFilename_or_aFilename (str or list): Single raster file path or list of tiled raster files
            load_data (bool): Whether to load raster data array (memory intensive)
            use_cache (bool): Whether to use/create cached instances
            iFlag_verbose_in (bool): Enable verbose logging

        Returns:
            sraster or list: Single sraster instance or list of sraster instances

        Note:
            - Supports both single files and tiled datasets (list of files)
            - Small files (< 100MB): Cache metadata for fast repeated access
            - Large files (> 100MB): Create transient instances to prevent memory issues
            - For tiled datasets: Efficiently handles metadata caching across all tiles
        """
        # Handle single file case
        if isinstance(sFilename_or_aFilename, str):
            return self._get_single_sraster(
                sFilename_or_aFilename, load_data, use_cache, iFlag_verbose_in
            )

        # Handle list of files (tiled dataset)
        elif isinstance(sFilename_or_aFilename, (list, tuple)):
            if iFlag_verbose_in:
                logger.info(
                    f"Processing tiled dataset with {len(sFilename_or_aFilename)} raster files"
                )

            sraster_list = []
            for idx, sFilename in enumerate(sFilename_or_aFilename):
                if iFlag_verbose_in:
                    logger.info(
                        f"  Processing tile {idx+1}/{len(sFilename_or_aFilename)}: {os.path.basename(sFilename)}"
                    )

                pRaster = self._get_single_sraster(
                    sFilename, load_data, use_cache, iFlag_verbose_in
                )
                sraster_list.append(pRaster)

            if iFlag_verbose_in:
                cached_count = sum(
                    1 for f in sFilename_or_aFilename if f in self._raster_cache
                )
                logger.info(
                    f"Tiled dataset processed: {cached_count}/{len(sFilename_or_aFilename)} tiles cached"
                )

            return sraster_list

        else:
            raise ValueError(
                f"sFilename_or_aFilename must be str or list, got {type(sFilename_or_aFilename)}"
            )

    def _get_single_sraster(
        self, sFilename, load_data=False, use_cache=True, iFlag_verbose_in=False
    ):
        """
        Get single sraster instance with intelligent caching.

        Internal method used by _get_sraster for individual file processing.
        """
        if not self._cache_enabled or not use_cache:
            # Cache disabled or explicitly requested not to cache
            if iFlag_verbose_in:
                logger.debug(
                    f"Creating transient sraster instance for: {os.path.basename(sFilename)}"
                )
            pRaster = sraster(sFilename)
            pRaster.read_metadata()
            return pRaster

        # Check if we already have this file cached
        if sFilename in self._raster_cache:
            if iFlag_verbose_in:
                logger.debug(
                    f"Using cached sraster metadata for: {os.path.basename(sFilename)}"
                )
            return self._raster_cache[sFilename]

        # File not cached, need to decide whether to cache it
        try:
            file_size = os.path.getsize(sFilename) if os.path.exists(sFilename) else 0
            should_cache = file_size < self._cache_size_threshold

            if iFlag_verbose_in:
                size_mb = file_size / (1024 * 1024)
                threshold_mb = self._cache_size_threshold / (1024 * 1024)
                logger.debug(
                    f"Raster {os.path.basename(sFilename)}: {size_mb:.1f}MB "
                    f"(threshold: {threshold_mb:.0f}MB, will_cache: {should_cache})"
                )

            pRaster = sraster(sFilename)
            pRaster.read_metadata()

            if should_cache and not load_data:
                # Cache metadata only for small files
                self._raster_cache[sFilename] = pRaster
                if iFlag_verbose_in:
                    logger.debug(f"Cached metadata for: {os.path.basename(sFilename)}")
            elif should_cache and load_data:
                if iFlag_verbose_in:
                    logger.debug(
                        f"Small file but data requested - not caching to avoid memory bloat: {os.path.basename(sFilename)}"
                    )
            else:
                if iFlag_verbose_in:
                    logger.debug(
                        f"Large file - using transient instance: {os.path.basename(sFilename)}"
                    )

            return pRaster

        except Exception as e:
            logger.warning(f"Error in _get_single_sraster for {sFilename}: {e}")
            # Fallback to simple creation
            pRaster = sraster(sFilename)
            pRaster.read_metadata()
            return pRaster

    def _clear_raster_cache(self):
        """Clear the raster cache to free memory."""
        if self._raster_cache:
            cache_count = len(self._raster_cache)
            self._raster_cache.clear()
            logger.info(f"Cleared raster cache ({cache_count} instances freed)")

    def _get_cache_info(self):
        """Get information about current cache state."""
        return {
            "cached_files": len(self._raster_cache),
            "cache_enabled": self._cache_enabled,
            "size_threshold_mb": self._cache_size_threshold / (1024 * 1024),
            "cached_filenames": [
                os.path.basename(f) for f in self._raster_cache.keys()
            ],
        }

    def get_sraster_with_data(self, sFilename_or_aFilename, iFlag_verbose_in=False):
        """
        Context manager for data-intensive sraster operations.

        Ensures proper cleanup of memory-intensive operations by creating
        transient instances that are automatically cleaned up after use.

        Args:
            sFilename_or_aFilename (str or list): Single file or list of raster files
            iFlag_verbose_in (bool): Enable verbose logging

        Returns:
            context manager: Yields sraster instance(s) with data loaded

        Example:
            with self.get_sraster_with_data('large_file.tif') as pRaster:
                pRaster.read_data()  # Load pixel data
                # Process data here
                # Automatic cleanup when exiting context
        """

        class SRasterDataContext:
            def __init__(self, uraster_instance, filename_or_list, verbose):
                self.uraster_instance = uraster_instance
                self.filename_or_list = filename_or_list
                self.verbose = verbose
                self.sraster_instances = None

            def __enter__(self):
                # Always create fresh instances for data operations, never cache data
                self.sraster_instances = self.uraster_instance._get_sraster(
                    self.filename_or_list,
                    load_data=False,  # Don't preload data, let user control when
                    use_cache=False,  # Never cache data-intensive instances
                    iFlag_verbose_in=self.verbose,
                )
                return self.sraster_instances

            def __exit__(self, exc_type, exc_val, exc_tb):
                # Cleanup: explicitly release any data arrays if they exist
                if isinstance(self.sraster_instances, list):
                    for pRaster in self.sraster_instances:
                        if hasattr(pRaster, "aData") and pRaster.aData is not None:
                            pRaster.aData = None
                        if (
                            hasattr(pRaster, "pDataset")
                            and pRaster.pDataset is not None
                        ):
                            pRaster.pDataset = None
                else:
                    pRaster = self.sraster_instances
                    if hasattr(pRaster, "aData") and pRaster.aData is not None:
                        pRaster.aData = None
                    if hasattr(pRaster, "pDataset") and pRaster.pDataset is not None:
                        pRaster.pDataset = None

                if self.verbose:
                    logger.debug("Released data arrays for memory cleanup")

        return SRasterDataContext(self, sFilename_or_aFilename, iFlag_verbose_in)

    def setup(self, iFlag_verbose_in=False):
        """
        Initialize and validate the uraster configuration.
        Checks raster files and mesh file for existence and validity.

        Args:
            iFlag_verbose_in (bool, optional): If True, print detailed progress messages.
                If False, only print error messages. Default is False.

        Returns:
            bool: True if setup successful, False otherwise
        """
        raster_check = self.check_raster_files(iFlag_verbose_in=iFlag_verbose_in)
        mesh_check = self.check_mesh_file(iFlag_verbose_in=iFlag_verbose_in)

        return raster_check is not None and mesh_check is not None

    def check_raster_files(
        self, aFilename_source_raster_in=None, iFlag_verbose_in=False
    ):
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
            iFlag_verbose_in (bool, optional): If True, print detailed progress messages.
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
            logger.error("No raster files provided for validation")
            return None

        if not isinstance(aFilename_source_raster, (list, tuple)):
            logger.error(
                f"Raster files must be provided as a list, got {type(aFilename_source_raster).__name__}"
            )
            return None

        if iFlag_verbose_in:
            logger.info(f"Validating {len(aFilename_source_raster)} raster file(s)...")

        # Phase 1: Check file existence and readability
        for idx, sFilename_raster_in in enumerate(aFilename_source_raster, 1):
            if not isinstance(sFilename_raster_in, str):
                logger.error(
                    f"Raster file path must be a string, got {type(sFilename_raster_in).__name__} at index {idx}"
                )
                return None

            if not sFilename_raster_in.strip():
                logger.error(f"Empty raster file path at index {idx}")
                return None

            if not os.path.exists(sFilename_raster_in):
                logger.error(f"Raster file does not exist: {sFilename_raster_in}")
                return None

            if not os.path.isfile(sFilename_raster_in):
                logger.error(f"Path is not a file: {sFilename_raster_in}")
                return None

            # Check file permissions
            if not os.access(sFilename_raster_in, os.R_OK):
                logger.error(f"Raster file is not readable: {sFilename_raster_in}")
                return None

            # Quick GDAL format validation
            try:
                pDataset_test = gdal.Open(sFilename_raster_in, gdal.GA_ReadOnly)
                if pDataset_test is None:
                    logger.error(f"GDAL cannot open raster file: {sFilename_raster_in}")
                    return None
                pDataset_test = None  # Close dataset
            except Exception as e:
                logger.error(
                    f"Error opening raster with GDAL: {sFilename_raster_in}: {e}"
                )
                return None

        if iFlag_verbose_in:
            logger.info("All raster files exist and are readable")

        # Phase 2: Process and convert rasters to WGS84
        aFilename_source_raster_out = []

        # Create WGS84 spatial reference for comparison
        pSpatialRef_wgs84 = None
        try:
            pSpatialRef_wgs84 = osr.SpatialReference()
            pSpatialRef_wgs84.ImportFromEPSG(4326)
            wkt_wgs84 = pSpatialRef_wgs84.ExportToWkt()
        except Exception as e:
            logger.error(f"Failed to create WGS84 spatial reference: {e}")
            return None
        finally:
            # Clean up spatial reference object
            if pSpatialRef_wgs84 is not None:
                pSpatialRef_wgs84 = None

        # Process each raster file
        aExtent = list()
        for idx, sFilename_raster_in in enumerate(aFilename_source_raster, 1):
            if iFlag_verbose_in:
                logger.info(
                    f"Processing raster {idx}/{len(aFilename_source_raster)}: {os.path.basename(sFilename_raster_in)}"
                )

            try:
                # Create sraster instance using intelligent caching and read metadata
                pRaster = self._get_sraster(
                    sFilename_raster_in,
                    load_data=False,
                    use_cache=True,
                    iFlag_verbose_in=iFlag_verbose_in,
                )

                # Validate critical metadata
                if pRaster.pSpatialRef_wkt is None:
                    logger.error(
                        f"Raster has no spatial reference: {sFilename_raster_in}"
                    )
                    return None

                if pRaster.nrow is None or pRaster.ncolumn is None:
                    logger.error(f"Invalid raster dimensions: {sFilename_raster_in}")
                    return None

                if pRaster.nrow <= 0 or pRaster.ncolumn <= 0:
                    logger.error(
                        f"Raster has invalid dimensions ({pRaster.nrow}x{pRaster.ncolumn}): {sFilename_raster_in}"
                    )
                    return None

                # Check if coordinate system matches WGS84
                if pRaster.pSpatialRef_wkt == wkt_wgs84:
                    if iFlag_verbose_in:
                        logger.info(f"  ✓ Already in WGS84 (EPSG:4326)")
                    aFilename_source_raster_out.append(sFilename_raster_in)
                    pExtent = pRaster.aExtent_wgs84
                    aExtent.append(pExtent)
                else:
                    # Convert to WGS84
                    if iFlag_verbose_in:
                        logger.info(
                            f'  → Converting to WGS84 from {pRaster.pSpatialRef.GetName() if pRaster.pSpatialRef else "unknown CRS"}'
                        )
                    try:
                        pRaster_wgs84 = pRaster.convert_to_wgs84()

                        if pRaster_wgs84 is None or not hasattr(
                            pRaster_wgs84, "sFilename"
                        ):
                            logger.error(
                                f"Conversion to WGS84 failed: {sFilename_raster_in}"
                            )
                            return None

                        if not os.path.exists(pRaster_wgs84.sFilename):
                            logger.error(
                                f"Converted WGS84 file not found: {pRaster_wgs84.sFilename}"
                            )
                            return None

                        if iFlag_verbose_in:
                            logger.info(f"  ✓ Converted to: {pRaster_wgs84.sFilename}")
                        aFilename_source_raster_out.append(pRaster_wgs84.sFilename)
                        pExtent = pRaster_wgs84.aExtent_wgs84
                        aExtent.append(pExtent)

                    except Exception as e:
                        logger.error(
                            f"Error during WGS84 conversion: {sFilename_raster_in}: {e}"
                        )
                        logger.error(f"Traceback: {traceback.format_exc()}")
                        return None

                # Log raster summary
                if iFlag_verbose_in:
                    logger.debug(
                        f"  - Dimensions: {pRaster.nrow} x {pRaster.ncolumn} pixels"
                    )
                    logger.debug(f"  - Data type: {pRaster.eType}")
                    if hasattr(pRaster, "dNoData"):
                        logger.debug(f"  - NoData value: {pRaster.dNoData}")

                pRaster.pSpatialRef = None  # Clean up spatial reference

            except AttributeError as e:
                logger.error(
                    f"Missing expected attribute in sraster: {sFilename_raster_in}: {e}"
                )
                logger.error(
                    f"Ensure sraster class has all required methods and attributes"
                )
                return None

            except Exception as e:
                logger.error(
                    f"Unexpected error processing raster {sFilename_raster_in}: {e}"
                )
                logger.error(f"Error type: {type(e).__name__}")
                logger.error(f"Traceback: {traceback.format_exc()}")
                return None

        # Final validation
        if len(aFilename_source_raster_out) != len(aFilename_source_raster):
            logger.error(
                f"Output count mismatch: expected {len(aFilename_source_raster)}, got {len(aFilename_source_raster_out)}"
            )
            return None

        # get overall extent
        if len(aExtent) > 0:
            xmin = min([extent[0] for extent in aExtent])
            xmax = max([extent[1] for extent in aExtent])
            ymin = min([extent[2] for extent in aExtent])
            ymax = max([extent[3] for extent in aExtent])
            self.aExtent_rasters = [xmin, xmax, ymin, ymax]
            if iFlag_verbose_in:
                logger.info(
                    f"Overall raster extent (WGS84): [{xmin}, {xmax}, {ymin}, {ymax}]"
                )

        if iFlag_verbose_in:
            logger.info(
                f"Successfully validated and prepared {len(aFilename_source_raster_out)} raster file(s)"
            )

        # Handle cache invalidation when raster list changes due to reprojection
        self._handle_raster_list_update(
            aFilename_source_raster, aFilename_source_raster_out, iFlag_verbose_in
        )

        self.aFilename_source_raster = aFilename_source_raster_out
        return aFilename_source_raster_out

    def _handle_raster_list_update(
        self, aFilename_original, aFilename_processed, iFlag_verbose_in=False
    ):
        """
        Handle cache invalidation when raster file list changes due to reprojection.

        When rasters are reprojected to WGS84, the output filenames differ from input filenames.
        This invalidates cached instances for the original files and requires cache management.

        Args:
            aFilename_original (list): Original input raster filenames
            aFilename_processed (list): Processed/reprojected raster filenames
            iFlag_verbose_in (bool): Enable verbose logging
        """
        if not self._cache_enabled or not self._raster_cache:
            return

        files_changed = aFilename_original != aFilename_processed

        if files_changed:
            if iFlag_verbose_in:
                logger.info(
                    "Raster file list changed due to reprojection - managing cache..."
                )

            # Count how many files were reprojected
            original_set = set(aFilename_original)
            processed_set = set(aFilename_processed)
            reprojected_count = len(original_set - processed_set)

            if reprojected_count > 0:
                if iFlag_verbose_in:
                    logger.info(f"  - {reprojected_count} file(s) were reprojected")

                # Remove cache entries for original files that were reprojected
                cache_keys_to_remove = []
                for original_file in aFilename_original:
                    if (
                        original_file not in processed_set
                        and original_file in self._raster_cache
                    ):
                        cache_keys_to_remove.append(original_file)

                for key in cache_keys_to_remove:
                    del self._raster_cache[key]
                    if iFlag_verbose_in:
                        logger.debug(
                            f"  - Removed cache entry for: {os.path.basename(key)}"
                        )

                if cache_keys_to_remove:
                    logger.info(
                        f"  - Cleared {len(cache_keys_to_remove)} stale cache entries"
                    )

                # Pre-cache metadata for new reprojected files if they're small enough
                if iFlag_verbose_in:
                    logger.info("  - Pre-caching metadata for reprojected files...")

                for processed_file in aFilename_processed:
                    if (
                        processed_file not in aFilename_original
                    ):  # This is a new reprojected file
                        try:
                            # Use _get_sraster to apply intelligent caching logic
                            self._get_sraster(
                                processed_file,
                                load_data=False,
                                use_cache=True,
                                iFlag_verbose_in=False,
                            )
                        except Exception as e:
                            if iFlag_verbose_in:
                                logger.warning(
                                    f"  - Could not pre-cache {os.path.basename(processed_file)}: {e}"
                                )

        else:
            if iFlag_verbose_in:
                logger.debug("No raster reprojection needed - cache remains valid")

    def check_mesh_file(self, iFlag_verbose_in=False):
        """
        Check if the source mesh file exists and build its topology.

        Args:
            iFlag_verbose_in (bool, optional): If True, print detailed progress messages.
                If False, only print error messages. Default is False.

        Returns:
            tuple or None: (vertices_lon, vertices_lat, connectivity) if successful, None otherwise
        """
        if not self.sFilename_source_mesh:
            logger.error("No source mesh filename provided")
            return None

        if not os.path.exists(self.sFilename_source_mesh):
            logger.error(
                f"Source mesh file does not exist: {self.sFilename_source_mesh}"
            )
            return None

        # need to add a function to check whether cellid field exists in the mesh file, if not, add it
        if self.sField_unique_id is None:
            logger.error(
                "No unique ID field specified for mesh cells (sField_unique_id)"
            )
            return None

        # check mesh geometry_validity, if there are invalid range, return None
        # if there are cell cross IDL, need to split the cell into two parts, but it will valid
        self.sFilename_source_mesh = self.check_mesh_geometry_validity(
            iFlag_verbose_in=iFlag_verbose_in
        )
        if not self.sFilename_source_mesh:
            return None

        # Setup and validate the mesh cell ID field
        if not self.setup_mesh_cellid():
            return None

        return self.rebuild_mesh_topology(iFlag_verbose_in=iFlag_verbose_in)

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
            ogr.wkbGeometryCollection: "wkbGeometryCollection",
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

    def check_mesh_geometry_validity(self, iFlag_verbose_in=False):

        return utility.check_mesh_quality(
            self.sFilename_source_mesh, iFlag_verbose_in=iFlag_verbose_in
        )

    def setup_mesh_cellid(self):
        """
        Setup and validate integer mesh cell ID field for GeoVista compatibility.

        Always ensures integer cell IDs by generating globally unique sequential integers.
        If field doesn't exist or is string type, creates intermediate file with new integer field.
        Original string field values are preserved when creating new integer IDs.

        Returns:
            bool: True if setup successful, False otherwise

        Note:
            - Always generates globally unique integer cell IDs (1, 2, 3, ...)
            - Preserves original field values when converting from string
            - Required for GeoVista visualization compatibility
        """
        if not self.sFilename_source_mesh or not self.sField_unique_id:
            logger.error("Missing mesh filename or unique ID field name")
            return False

        try:
            # Open mesh file to check field existence
            pDataset_mesh = ogr.Open(self.sFilename_source_mesh, 0)  # Read-only
            if pDataset_mesh is None:
                logger.error(f"Cannot open mesh file: {self.sFilename_source_mesh}")
                return False

            pLayer_mesh = pDataset_mesh.GetLayer()
            if pLayer_mesh is None:
                logger.error("Cannot access layer in mesh file")
                pDataset_mesh = None
                return False

            # Check if specified field exists
            pLayer_defn = pLayer_mesh.GetLayerDefn()
            iFlag_field_exists = False

            for i in range(pLayer_defn.GetFieldCount()):
                pField_defn = pLayer_defn.GetFieldDefn(i)
                if pField_defn.GetName().lower() == self.sField_unique_id.lower():
                    iFlag_field_exists = True
                    iField_type = pField_defn.GetType()

                    # Always require integer type for GeoVista compatibility
                    if iField_type == ogr.OFTInteger:
                        self.iField_unique_type = ogr.OFTInteger
                        logger.info(
                            f"Found integer field '{self.sField_unique_id}' - ready for use"
                        )
                        pDataset_mesh = None
                        return True
                    else:
                        # Field exists but wrong type - regenerate with integer IDs
                        logger.info(
                            f"Found field '{self.sField_unique_id}' but not integer type - will regenerate with sequential integer IDs"
                        )
                        iFlag_field_exists = False  # Treat as needing regeneration
                    break

            if not iFlag_field_exists:
                # Field doesn't exist or needs integer regeneration
                logger.info(
                    f"Creating intermediate file with sequential integer field '{self.sField_unique_id}'..."
                )

                # Create intermediate filename
                base_name = os.path.splitext(self.sFilename_source_mesh)[0]
                extension = os.path.splitext(self.sFilename_source_mesh)[1]
                sFilename_intermediate = (
                    f"{base_name}_with_{self.sField_unique_id}{extension}"
                )

                # Copy original file and add the field
                success = self._create_intermediate_mesh_file(
                    self.sFilename_source_mesh,
                    sFilename_intermediate,
                    self.sField_unique_id,
                    pDataset_mesh,
                )

                pDataset_mesh = None

                if success:
                    # Update the source mesh filename to use intermediate file
                    self.sFilename_source_mesh = sFilename_intermediate
                    logger.info(
                        f"Successfully created intermediate file: {sFilename_intermediate}"
                    )
                    logger.info(
                        f"Updated source mesh path to: {self.sFilename_source_mesh}"
                    )
                    return True
                else:
                    logger.error("Failed to create intermediate mesh file")
                    return False

        except Exception as e:
            logger.error(f"Error setting up mesh cell ID field: {e}")
            return False

    def _create_intermediate_mesh_file(
        self, sFilename_source, sFilename_target, sField_name, pDataset_source
    ):
        """
        Create intermediate mesh file with the specified ID field added.

        Args:
            sFilename_source (str): Path to source mesh file
            sFilename_target (str): Path to target intermediate file
            sField_name (str): Name of the ID field to add
            pDataset_source: Open source dataset (will be closed and reopened)

        Returns:
            bool: True if successful, False otherwise
        """
        try:
            # Close source dataset if open
            if pDataset_source is not None:
                pDataset_source = None

            # Reopen source dataset
            pDataset_source = ogr.Open(sFilename_source, 0)  # Read-only
            if pDataset_source is None:
                return False

            pLayer_source = pDataset_source.GetLayer()
            if pLayer_source is None:
                pDataset_source = None
                return False

            # Get source layer info
            pSpatialRef = pLayer_source.GetSpatialRef()
            geometry_type = pLayer_source.GetGeomType()

            # Create target dataset using the same driver as source
            pDriver = get_vector_driver_from_filename(sFilename_target)
            if os.path.exists(sFilename_target):
                pDriver.DeleteDataSource(sFilename_target)

            pDataset_target = pDriver.CreateDataSource(sFilename_target)
            if pDataset_target is None:
                pDataset_source = None
                return False

            # Create target layer with same properties
            pLayer_target = pDataset_target.CreateLayer(
                "mesh", pSpatialRef, geometry_type
            )
            if pLayer_target is None:
                pDataset_source = None
                pDataset_target = None
                return False

            # Copy all existing fields
            pLayerDefn_source = pLayer_source.GetLayerDefn()
            for i in range(pLayerDefn_source.GetFieldCount()):
                pFieldDefn = pLayerDefn_source.GetFieldDefn(i)
                pLayer_target.CreateField(pFieldDefn)

            # Check if we need to preserve existing field values
            pLayerDefn_source = pLayer_source.GetLayerDefn()
            iSource_field_index = pLayerDefn_source.GetFieldIndex(sField_name)
            bHas_existing_field = iSource_field_index >= 0

            # If field exists, add "oldid" field to preserve original values
            if bHas_existing_field:
                pField_defn_source = pLayerDefn_source.GetFieldDefn(iSource_field_index)
                iOriginal_field_type = pField_defn_source.GetType()
                pField_oldid = ogr.FieldDefn("oldid", iOriginal_field_type)
                pLayer_target.CreateField(pField_oldid)

            # Add the new integer ID field (will be stored globally)
            pField_id = ogr.FieldDefn(sField_name, ogr.OFTInteger)
            pLayer_target.CreateField(pField_id)
            # Set field type globally for newly created field
            self.iField_unique_type = ogr.OFTInteger

            # Copy features and add sequential IDs
            pLayer_source.ResetReading()
            feature_id = 1

            for pFeature_source in pLayer_source:
                # Create new feature
                pFeature_target = ogr.Feature(pLayer_target.GetLayerDefn())

                # Copy geometry
                pGeometry = pFeature_source.GetGeometryRef()
                if pGeometry is not None:
                    pFeature_target.SetGeometry(pGeometry.Clone())

                # Copy all existing field values and preserve original ID if exists
                for i in range(pLayerDefn_source.GetFieldCount()):
                    pFieldDefn = pLayerDefn_source.GetFieldDefn(i)
                    sField_name_existing = pFieldDefn.GetName()
                    field_value = pFeature_source.GetField(sField_name_existing)

                    # If this is the field being replaced, save to "oldid"
                    if (
                        sField_name_existing.lower() == sField_name.lower()
                        and bHas_existing_field
                    ):
                        pFeature_target.SetField("oldid", field_value)
                    elif sField_name_existing.lower() != sField_name.lower():
                        # Copy other fields normally (skip field being replaced)
                        pFeature_target.SetField(sField_name_existing, field_value)

                # Set the new globally unique integer ID field
                pFeature_target.SetField(sField_name, feature_id)

                # Add feature to target layer
                pLayer_target.CreateFeature(pFeature_target)

                pFeature_target = None
                feature_id += 1

            # Cleanup
            pDataset_source = None
            pDataset_target = None

            logger.info(f"Added field '{sField_name}' to {feature_id} features")
            return True

        except Exception as e:
            logger.error(f"Error creating intermediate file: {e}")
            logger.error(f"Traceback: {traceback.format_exc()}")
            return False

    def rebuild_mesh_topology(self, iFlag_verbose_in=False):
        """
        Rebuild mesh topology from source mesh file by extracting vertices,
        connectivity, and centroids for unstructured mesh processing.

        This method uses the enhanced standalone rebuild_mesh_topology function
        and updates all instance attributes with the comprehensive mesh information.

        Args:
            iFlag_verbose_in (bool, optional): If True, print detailed progress messages.
                If False, only print error messages. Default is False.

        Returns:
            tuple: (vertices_longitude, vertices_latitude, connectivity) or None on failure
        """
        # Use the enhanced standalone function from _visual module
        mesh_info = utility.rebuild_mesh_topology(
            self.sFilename_source_mesh,
            iFlag_verbose_in=iFlag_verbose_in,
            sField_unique_id=self.sField_unique_id,
        )

        if mesh_info is None:
            return None

        # Update all instance attributes with comprehensive mesh information
        self.aVertex_longititude = mesh_info["vertices_longitude"]
        self.aVertex_latitude = mesh_info["vertices_latitude"]
        self.aConnectivity = mesh_info["connectivity"]
        self.aCenter_longititude = mesh_info["cell_centroids_longitude"]
        self.aCenter_latitude = mesh_info["cell_centroids_latitude"]
        self.aCellID = mesh_info[
            "cell_ids"
        ]  # can be both int or str depending on the field type
        self.dArea_min = mesh_info["area_min"]
        self.dArea_max = mesh_info["area_max"]
        self.dArea_mean = mesh_info["area_mean"]
        self.nVertex_max = mesh_info["max_vertices_per_cell"]
        self.nPolygon = mesh_info["num_polygns"]  # Update cell count
        self.nCell_source = mesh_info["num_cells"]

        # Return the traditional tuple format for backward compatibility
        return (
            mesh_info["vertices_longitude"],
            mesh_info["vertices_latitude"],
            mesh_info["connectivity"],
        )

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

                print("\n" + "=" * 60)
                print("GPU/GeoVista Information:")
                print("=" * 60)
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
                f"Output file size: {os.path.getsize(sFilename_output) / (1024*1024):.2f} MB"
            )
        else:
            logger.warning("No output file information available")

    def print_raster_info(self):
        """
        Print detailed information about all input raster files.
        """
        print("\n" + "=" * 60)
        print(
            f"Input Raster Information ({len(self.aFilename_source_raster)} file(s)):"
        )
        print("=" * 60)

        for idx, sFilename in enumerate(self.aFilename_source_raster, 1):
            print(f"\n[{idx}] {sFilename}")
            try:
                # Use cached sraster instance for metadata access
                pRaster = self._get_sraster(
                    sFilename, load_data=False, use_cache=True, iFlag_verbose_in=False
                )
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

        print("\n" + "=" * 60)
        print("Mesh Topology Information:")
        print("=" * 60)
        print(f"Number of mesh cells: {len(self.aCenter_longititude)}")
        print(
            f"Cell longitude range: {self.aCenter_longititude.min():.3f} to {self.aCenter_longititude.max():.3f}"
        )
        print(
            f"Cell latitude range: {self.aCenter_latitude.min():.3f} to {self.aCenter_latitude.max():.3f}"
        )
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
                f"  Area range ratio: {self.dArea_max/self.dArea_min:.2f}x"
                if self.dArea_min > 0
                else "  Area range ratio: N/A"
            )

        print("=" * 60)

    def run_remap(
        self,
        sFilename_target_mesh_out=None,
        sFilename_source_mesh_in=None,
        aFilename_source_raster_in=None,
        iFlag_stat_in=True,
        iFlag_weighted_average_in=False,
        iFlag_remap_method_in=1,
        iFlag_save_clipped_raster_in=0,
        sFolder_raster_out_in=None,
        iFlag_discrete_in=False,
        iFlag_verbose_in=False,
    ):
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

        # check stat and discrete compatibility
        if iFlag_discrete_in:
            # for discrete, only remap method 1 (nearest neighbor) is allowed
            # can we apply statistics for discrete?
            self.iFlag_discrete = True
            iFlag_stat_in = False
            if iFlag_remap_method_in != 1:
                logger.error(
                    "For discrete remap, only remap method 1 (nearest neighbor) is allowed."
                )
                return None
        else:
            # for continuous, all remap methods are allowed
            pass

        # the model should suport weighted average and discrete remap
        if iFlag_weighted_average_in:
            # call the polygon calculation with weighted average
            sFilename_raster = aFilename_source_raster[0]
            # Use cached sraster instance for metadata access
            pRaster = self._get_sraster(
                sFilename_raster,
                load_data=False,
                use_cache=True,
                iFlag_verbose_in=iFlag_verbose_in,
            )
            sFilename_raster_mesh = pRaster.create_raster_mesh()
            return intersect.run_remap(
                sFilename_target_mesh,
                sFilename_source_mesh,
                sFilename_raster,
                sFilename_raster_mesh,
                iFlag_save_clipped_raster_in=iFlag_save_clipped_raster_in,
                sFolder_raster_out_in=sFolder_raster_out_in,
                iFlag_discrete_in=iFlag_discrete_in,
                iFlag_verbose_in=iFlag_verbose_in,
            )

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
                iFlag_verbose_in=iFlag_verbose_in,
                sField_unique_id=self.sField_unique_id,
            )

    def visualize_raster(self, sFilename_out=None, iFlag_verbose_in=False):
        """
        Visualize source raster data using GeoVista.

        Note:
            Not yet implemented. Placeholder for future raster visualization.
        """
        return _visual.visualize_raster(
            self, sFilename_out=sFilename_out, iFlag_verbose_in=iFlag_verbose_in
        )

    def visualize_source_mesh(self, sFilename_out=None, **kwargs):
        """
        Visualize the source mesh topology using GeoVista 3D globe rendering.

        Creates an interactive or saved 3D visualization of the unstructured mesh
        with proper geographic context including coastlines and coordinate grid.

        Args:
            sFilename_out (str, optional): Output screenshot file path.
                If None, displays interactive viewer. Supports formats: .png, .jpg, .svg

        Keyword Arguments (all optional with sensible defaults):
            dLongitude_focus_in (float): Camera focal point longitude in degrees.
                Valid range: -180 to 180. Default is 0.0 (prime meridian).
            dLatitude_focus_in (float): Camera focal point latitude in degrees.
                Valid range: -90 to 90. Default is 0.0 (equator).
            dImage_scale_in (float): Image scaling factor. Default is 1.0.
            dZoom_factor (float): Camera zoom level. Higher values zoom in. Default is 0.7.
            window_size_in (tuple): Window size as (width, height). Default is (800, 600).
            iFlag_show_coastlines (bool): Show coastline overlay. Default is True.
            iFlag_show_graticule (bool): Show coordinate grid with labels. Default is True.
            iFlag_wireframe_only (bool): Show wireframe only. Default is True.
            iFlag_verbose_in (bool): If True, print detailed progress messages. Default is False.

        Returns:
            bool: True if visualization successful, False otherwise

        Example:
            # Simple call with defaults
            obj.visualize_source_mesh()

            # Override specific parameters
            obj.visualize_source_mesh('output.png', dZoom_factor=0.5, iFlag_show_coastlines=False)

            # Customize multiple parameters
            obj.visualize_source_mesh(sFilename_out='map.png',
                                      dLongitude_focus_in=120.0,
                                      dLatitude_focus_in=30.0)

        Note:
            - Requires 'geovista' package: pip install geovista
            - Interactive mode requires display environment
            - Mesh topology must be built before visualization (call rebuild_mesh_topology first)
        """
        # Set defaults for all parameters
        defaults = {
            "dLongitude_focus_in": 0.0,
            "dLatitude_focus_in": 0.0,
            "dImage_scale_in": 1.0,
            "dZoom_factor": 0.7,
            "window_size_in": (800, 600),
            "iFlag_show_coastlines": True,
            "iFlag_show_graticule": True,
            "iFlag_wireframe_only": True,
            "iFlag_verbose_in": False,
        }

        # Merge defaults with provided kwargs
        merged_params = {**defaults, **kwargs}

        return _visual.visualize_source_mesh(self, sFilename_out, **merged_params)

    def visualize_target_mesh(
        self, sVariable_in=None, sUnit_in=None, sFilename_out=None, **kwargs
    ):
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

        Keyword Arguments (all optional with sensible defaults):
            dLongitude_focus_in (float): Camera focal point longitude in degrees.
                Valid range: -180 to 180. Default is 0.0. For animations, this is the starting longitude.
            dLatitude_focus_in (float): Camera focal point latitude in degrees.
                Valid range: -90 to 90. Default is 0.0.
            dImage_scale_in (float): Image scaling factor. Default is 1.0.
            dZoom_factor (float): Camera zoom level. Higher values zoom in. Default is 0.7.
            window_size_in (tuple): Window size as (width, height). Default is (800, 600).
            iFlag_show_coastlines (bool): Show coastline overlay. Default is True.
            iFlag_show_graticule (bool): Show coordinate grid with labels. Default is True.
            sColormap (str): Matplotlib colormap name. Default is 'viridis'.
                Examples: 'plasma', 'coolwarm', 'jet', 'RdYlBu'
            iFlag_create_animation (bool): Create rotating animation. Default is False.
            iAnimation_frames (int): Number of frames for 360° rotation. Default is 360.
            dAnimation_speed (float): Animation speed in degrees per frame. Default is 1.0.
            sAnimation_format (str): Animation output format. Default is 'mp4'.
                Supports: 'mp4', 'gif', 'avi'
            iFlag_verbose_in (bool): If True, print detailed progress messages. Default is False.

        Returns:
            bool: True if visualization successful, False otherwise

        Example:
            # Simple visualization
            obj.visualize_target_mesh('mean')

            # Save to file with custom colormap
            obj.visualize_target_mesh('mean', sFilename_out='output.png', sColormap='coolwarm')

            # Create animation
            obj.visualize_target_mesh('mean', sFilename_out='animation.mp4',
                                      iFlag_create_animation=True, iAnimation_frames=120)

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
        # Set defaults for all parameters
        defaults = {
            "dLongitude_focus_in": 0.0,
            "dLatitude_focus_in": 0.0,
            "dImage_scale_in": 1.0,
            "dZoom_factor": 0.7,
            "window_size_in": (800, 600),
            "iFlag_show_coastlines": True,
            "iFlag_show_graticule": False,
            "sColormap": "viridis",
            "iFlag_create_animation": False,
            "iAnimation_frames": 360,
            "dAnimation_speed": 1.0,
            "sAnimation_format": "mp4",
            "iFlag_verbose_in": False,
        }

        # Merge defaults with provided kwargs
        merged_params = {**defaults, **kwargs}
        if sVariable_in is None:
            if self.iFlag_discrete:
                sVariable_in = "mode"
            else:
                sVariable_in = "mean"

        return _visual.visualize_target_mesh(
            self, sVariable_in, sUnit_in, sFilename_out, **merged_params
        )

    def _create_rotation_animation(
        self,
        plotter,
        sFilename_out,
        dLongitude_start,
        dLatitude_focus,
        iAnimation_frames,
        dAnimation_speed,
        sAnimation_format,
        iFlag_verbose_in=False,
    ):
        """
        Create a rotating animation of the 3D globe visualization.

        This method delegates to the _visual module for implementation.
        """
        return _visual._create_rotation_animation(
            self,
            plotter,
            sFilename_out,
            dLongitude_start,
            dLatitude_focus,
            iAnimation_frames,
            dAnimation_speed,
            sAnimation_format,
            iFlag_verbose_in,
        )

    def cleanup(self):
        """
        Comprehensive cleanup method to release spatial reference objects,
        cached raster instances, and other resources.
        """
        try:
            # Clean up raster cache
            if hasattr(self, "_raster_cache") and self._raster_cache:
                cache_count = len(self._raster_cache)
                self._raster_cache.clear()
                logger.debug(f"Cleared raster cache ({cache_count} instances)")

            # Clean up spatial reference objects
            if hasattr(self, "pSpatialRef") and self.pSpatialRef is not None:
                self.pSpatialRef = None
                logger.debug("Spatial reference object cleaned up successfully")

            # Reset cache configuration
            if hasattr(self, "_cache_enabled"):
                self._cache_enabled = True  # Reset to default state

        except Exception as e:
            logger.warning(f"Error during cleanup: {e}")

    def set_cache_enabled(self, enabled=True):
        """
        Enable or disable raster caching.

        Args:
            enabled (bool): Whether to enable caching. Default is True.

        Note:
            Disabling cache clears existing cached instances.
        """
        if not enabled and hasattr(self, "_raster_cache"):
            self._clear_raster_cache()
        self._cache_enabled = enabled
        logger.info(f"Raster caching {'enabled' if enabled else 'disabled'}")

    def set_cache_threshold(self, threshold_mb=100):
        """
        Set the file size threshold for caching decisions.

        Args:
            threshold_mb (int): File size threshold in megabytes. Default is 100MB.
        """
        self._cache_size_threshold = threshold_mb * 1024 * 1024
        logger.info(f"Cache size threshold set to {threshold_mb}MB")
