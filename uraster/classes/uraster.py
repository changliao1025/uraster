# Define a class named 'uraster'
import os
import logging
import traceback
import signal
import sys
import time
import numpy as np
from osgeo import gdal, ogr, osr
from multiprocessing import Pool, cpu_count
from concurrent.futures import ThreadPoolExecutor, as_completed
from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
from pyearth.gis.geometry.split_polygon_cross_idl import split_polygon_cross_idl
from .sraster import sraster

# Set up logging for crash detection
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def signal_handler(signum, frame):
    """Handle system signals (like SIGSEGV) for crash detection"""
    logger.error(f"Received signal {signum}. Process may have crashed.")
    logger.error(f"Current frame: {frame}")
    traceback.print_stack(frame)
    sys.exit(1)

# Register signal handlers for common crash signals
signal.signal(signal.SIGSEGV, signal_handler)  # Segmentation fault
signal.signal(signal.SIGABRT, signal_handler)  # Abort
signal.signal(signal.SIGFPE, signal_handler)   # Floating point exception

class uraster:

    def __init__(self, aConfig = dict()):
        self.iFlag_global = None
        self.dResolution_raster = None
        self.dResolution_uraster = None
        self.sFilename_source_mesh = None
        self.sFilename_target_mesh = None
        self.aFilename_source_raster = list()

        if "sFilename_source_mesh" in aConfig:
            self.sFilename_source_mesh = aConfig['sFilename_source_mesh']

        if "sFilename_target_mesh" in aConfig:
            self.sFilename_target_mesh = aConfig['sFilename_target_mesh']

        if "aFilename_source_raster" in aConfig:
            self.aFilename_source_raster = aConfig['aFilename_source_raster']

    def _get_geometry_type_name(self, geometry_type):
        """
        Convert OGR geometry type integer to readable string
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

    def _get_vector_driver_from_filename(self, filename):
        """
        Determine the appropriate OGR driver based on file extension
        :param filename: Output vector filename
        :return: OGR driver name
        """
        extension = os.path.splitext(filename)[1].lower()

        driver_mapping = {
            '.shp': 'ESRI Shapefile',
            '.geojson': 'GeoJSON',
            '.json': 'GeoJSON',
            '.gpkg': 'GPKG',
            '.kml': 'KML',
            '.gml': 'GML',
            '.sqlite': 'SQLite',
            '.csv': 'CSV',
            '.parquet': 'Parquet'
        }

        # Default to Shapefile if extension not recognized
        return driver_mapping.get(extension, 'ESRI Shapefile')

    def check_raster_files(self, aFilename_source_raster_in=None):
        """
        Check if the input raster files exist
        :param aFilename_source_raster_in: a list of raster files
        :return: True if all files exist, False otherwise
        """

        if aFilename_source_raster_in is None:
            aFilename_source_raster = self.aFilename_source_raster
        else:
            aFilename_source_raster = aFilename_source_raster_in

        for sFilename_raster_in in aFilename_source_raster:
            if os.path.exists(sFilename_raster_in):
                pass
            else:
                print('The raster file does not exist:', sFilename_raster_in)
                return False

        aFilename_source_raster_out = list()
        #create a WGS84 spatial reference in the WKT format for comparison
        pSpatialRef_wgs84 = osr.SpatialReference()
        pSpatialRef_wgs84.ImportFromEPSG(4326)
        wkt_wgs84 = pSpatialRef_wgs84.ExportToWkt()
        #uset the sraster class the check the raster
        for sFilename_raster_in in aFilename_source_raster:
            pRaster = sraster(sFilename_raster_in)
            pRaster.read_metadata()
            if pRaster.pSpatialRef_wkt == wkt_wgs84:
                print('The raster file is in WGS84 geographic coordinate system:', sFilename_raster_in)
                aFilename_source_raster_out.append(sFilename_raster_in)

            else:
                #need conversion
                pRaster_wgs84 = pRaster.convert_to_wgs84()
                aFilename_source_raster_out.append(pRaster_wgs84.sFilename)


        return aFilename_source_raster_out

    def print_raster_info(self):
        """
        Print the raster information
        :return: None
        """

        print('The input raster files are:')
        for sFilename in self.aFilename_source_raster:
            print(sFilename)
            pRaster = sraster(sFilename)
            pRaster.read_metadata()
            pRaster.print_info()

        return

    def remap_raster_to_uraster(self, sFilename_vector_out,
                    sFilename_target_mesh_in = None,
                                aFilename_source_raster_in = None,
                  iFlag_stat_in = 1,
                  iFlag_save_clipped_raster_in=0,
                  sFolder_raster_out_in = None,
                  sFormat_in='GTiff'):
        """
        Clip a raster by a vector mesh file
        :param aFilename_source_raster_in: a list of raster files
        :param    sFilename_target_mesh_in : input polygon filename
        :param sFilename_raster_out: output raster filename
        :param sFormat: output format
        :return: None
        """

        if aFilename_source_raster_in is None:
            aFilename_source_raster = self.aFilename_source_raster
        else:
            aFilename_source_raster = aFilename_source_raster_in

        if sFilename_target_mesh_in is None:
            sFilename_target_mesh = self.sFilename_target_mesh
        else:
            sFilename_target_mesh = sFilename_target_mesh_in

        #check input files
        for sFilename_raster in aFilename_source_raster:
            if os.path.exists(sFilename_raster):
                pass
            else:
                print('The raster file does not exist!')
                return

        if os.path.exists(sFilename_target_mesh ):
            pass
        else:
            print('The vector mesh file does not exist!')
            return

        # Determine output vector format from filename extension
        sVectorDriverName = self._get_vector_driver_from_filename(sFilename_vector_out)
        pDriver_vector = ogr.GetDriverByName(sVectorDriverName)

        #check the input raster data format and decide gdal driver
        if sFormat_in is not None:
            sDriverName = sFormat_in
        else:
            sDriverName = 'GTiff'

        if os.path.exists(sFilename_vector_out):
            #remove the file using the vector driver
            pDriver_vector.DeleteDataSource(sFilename_vector_out)

        pDriver = gdal.GetDriverByName(sDriverName)

        #get the raster file extension
        sFilename_raster = aFilename_source_raster[0]  #just use the first raster to get the extension
        sExtension = os.path.splitext(sFilename_raster)[1]
        sName = os.path.basename(sFilename_raster)
        sRasterName_no_extension = os.path.splitext(sName)[0]

        #use the sraster class the check the raster
        dX_left = -180.0
        dX_right = 180.0
        dY_top = 90.0
        dY_bot = -90.0
        dPixelWidth = None
        pPixelHeight = None
        #narrow the range to speed up the processing
        #also get the highest resolution
        for sFilename_raster in aFilename_source_raster:
            sRaster = sraster(sFilename_raster)
            sRaster.read_metadata()
            eType = sRaster.eType
            dX_left = min(dX_left, sRaster.dLongitude_left)
            dX_right = max(dX_right, sRaster.dLongitude_right)
            dY_top = max(dY_top, sRaster.dLatitude_top)
            dY_bot = min(dY_bot, sRaster.dLatitude_bottom)
            dMissing_value = sRaster.dNoData
            if dPixelWidth is None or sRaster.pTransform[1] < dPixelWidth:
                dPixelWidth = sRaster.pTransform[1]

            if pPixelHeight is None or sRaster.pTransform[5] > pPixelHeight: #is pPixelHeight is negative, then it should be less than
                pPixelHeight = sRaster.pTransform[5]


        #get the spatial reference of the mesh vector file
        pDataset_mesh = ogr.Open( sFilename_target_mesh )
        pLayer_mesh = pDataset_mesh.GetLayer(0)
        nFeature = pLayer_mesh.GetFeatureCount()
        pSpatialRef_target = pLayer_mesh.GetSpatialRef()

        pSpatialRef_target_wkt = pSpatialRef_target.ExportToWkt()

        #check whether the polygon has only one or more features
        if nFeature > 1:
            pass
        else:
            print('The polygon file has only one polygons!')
            return

        options = ['COMPRESS=DEFLATE', 'PREDICTOR=2']
        sResampleAlg= 'near'
        sDriverName = 'MEM'

        #create a polygon feature to save the output
        pDataset_out = pDriver_vector.CreateDataSource(sFilename_vector_out)
        pLayer_out = pDataset_out.CreateLayer('cell', pSpatialRef_target, ogr.wkbPolygon)
        pLayer_defn_out = pLayer_out.GetLayerDefn()
        pFeature_out = ogr.Feature(pLayer_defn_out)

        #add id, area and mean, min, max, std of the raster
        pLayer_out.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
        #define a field
        pField = ogr.FieldDefn('area', ogr.OFTReal)
        pField.SetWidth(32)
        pField.SetPrecision(2)
        pLayer_out.CreateField(pField)

        #in the future, we will also copy other attributes from the input geojson file

        if iFlag_stat_in == 1:
            pLayer_out.CreateField(ogr.FieldDefn('mean', ogr.OFTReal))
            pLayer_out.CreateField(ogr.FieldDefn('min', ogr.OFTReal))
            pLayer_out.CreateField(ogr.FieldDefn('max', ogr.OFTReal))
            pLayer_out.CreateField(ogr.FieldDefn('std', ogr.OFTReal))
        else:
            #treat the raster as categorical?
            pass

        # Pre-compute GDAL options to avoid repeated object creation
        gdal_warp_options_base = {
            'cropToCutline': True,
            'xRes': dPixelWidth,
            'yRes': abs(pPixelHeight),
            'dstSRS': pSpatialRef_target,
            'format': 'MEM',
            'resampleAlg': sResampleAlg,
            'srcSRS': 'EPSG:4326',  # Explicitly set source CRS
        }

        # Enable GDAL multi-threading
        gdal.SetConfigOption('GDAL_NUM_THREADS', str(cpu_count()))

        # Batch processing variables
        i = 1

        # Pre-fetch all features for potential parallel processing
        aFeatures = []
        pLayer_mesh.ResetReading()
        pFeature_mesh = pLayer_mesh.GetNextFeature()
        while pFeature_mesh is not None:
            #what if a feature cross the international date line?
            #reuse the code from pyflowline to split the polygon if it cross the date line
            pPolygon = pFeature_mesh.GetGeometryRef()
            aCoord = get_geometry_coordinates(pPolygon)
            dLon_min = np.min(aCoord[:,0])
            dLon_max = np.max(aCoord[:,0])
            if dLon_max - dLon_min > 100:  #cross the international date line
                if not pPolygon.IsValid():
                    print('Invalid polygon geometry detected, something broke during the PyFlowline simulation.')
                    pass
                # Create multipolygon to handle IDL crossing
                print('Feature ', i, ' cross the international date line, splitting it into multiple parts.')
                pMultipolygon = ogr.Geometry(ogr.wkbMultiPolygon)
                aCoord_gcs_split = split_polygon_cross_idl(aCoord) #careful
                for aCoord_gcs in aCoord_gcs_split:
                    #create a polygon (not just ring) and add it to the multipolygon
                    ring = ogr.Geometry(ogr.wkbLinearRing)
                    for iCoord in range(aCoord_gcs.shape[0]):
                        ring.AddPoint(aCoord_gcs[iCoord, 0], aCoord_gcs[iCoord, 1])

                    ring.CloseRings()
                    # Create polygon from ring
                    polygon_part = ogr.Geometry(ogr.wkbPolygon)
                    polygon_part.AddGeometry(ring)
                    pMultipolygon.AddGeometry(polygon_part)  # Add polygon, not ring

                #create a feature from the multipolygon and add it to the list
                pFeature_new = ogr.Feature(pLayer_mesh.GetLayerDefn())
                pFeature_new.SetGeometry(pMultipolygon)
                aFeatures.append(pFeature_new)
                pFeature_mesh = pLayer_mesh.GetNextFeature()
            else:
                aFeatures.append(pFeature_mesh.Clone())
                pFeature_mesh = pLayer_mesh.GetNextFeature()

            i += 1

        print(f"Processing {len(aFeatures)} features...")
        start_time = time.time()

        # Add crash tracking variables
        failed_features = []
        successful_features = 0

        #reset i
        i = 1

        # Process features
        for idx, pFeature_mesh in enumerate(aFeatures):
            sClip = f"{i:08d}"  # f-string is faster than format
            pPolygon = pFeature_mesh.GetGeometryRef()
            # Fast envelope check first
            minX, maxX, minY, maxY = pPolygon.GetEnvelope()
            if (minX > dX_right or maxX < dX_left or
                minY > dY_top or maxY < dY_bot or
                not pPolygon or pPolygon.IsEmpty() or not pPolygon.IsValid()):
                i += 1
                continue

            # Process valid polygons with comprehensive error handling
            feature_start_time = time.time()
            try:
                # Monitor memory usage (optional - requires psutil)
                try:
                    import psutil
                    process = psutil.Process()
                    memory_percent = process.memory_percent()
                    if memory_percent > 80:  # If using more than 80% of available memory
                        logger.warning(f"High memory usage: {memory_percent:.1f}% at feature {i}")
                except ImportError:
                    pass  # psutil not available

                # Handle polygon vs multipolygon differently
                geometry_type = pPolygon.GetGeometryType()
                sGeometry_type = pPolygon.GetGeometryName()
                geometry_type_name = self._get_geometry_type_name(geometry_type)
                logger.debug(f"Feature {i} geometry type: {geometry_type} ({geometry_type_name})")

                if sGeometry_type == "POLYGON":
                    # Simple polygon - process normally
                    pDataset_clip_warped, aData_clip, newGeoTransform = self._process_single_polygon(
                        pPolygon, aFilename_source_raster, gdal_warp_options_base, i)

                    if pDataset_clip_warped is None:
                        error_msg = f"Failed to process single polygon for feature {i}"
                        logger.error(error_msg)
                        failed_features.append({'feature_id': i, 'error': error_msg, 'envelope': (minX, maxX, minY, maxY)})
                        i += 1
                        continue

                elif sGeometry_type == "MULTIPOLYGON":
                    # Multipolygon (IDL crossing) - process each part separately and merge
                    logger.info(f"Processing IDL-crossing multipolygon for feature {i}")
                    merged_data, merged_transform = self._process_multipolygon_idl(
                        pPolygon, aFilename_source_raster, gdal_warp_options_base, i, dMissing_value)

                    if merged_data is None:
                        error_msg = f"Failed to process IDL-crossing multipolygon for feature {i}"
                        logger.error(error_msg)
                        failed_features.append({'feature_id': i, 'error': error_msg, 'envelope': (minX, maxX, minY, maxY)})
                        i += 1
                        continue

                    aData_clip = merged_data
                    newGeoTransform = merged_transform
                    pDataset_clip_warped = True  # Flag to indicate successful processing

                else:
                    error_msg = f"Unsupported geometry type for feature {i}: {geometry_type} ({geometry_type_name})"
                    logger.error(error_msg)
                    failed_features.append({'feature_id': i, 'error': error_msg, 'envelope': (minX, maxX, minY, maxY)})
                    i += 1
                    continue

                logger.debug(f"Successfully processed GDAL Warp for feature {i} in {time.time() - feature_start_time:.2f}s")
                successful_features += 1

            except TimeoutError as e:
                logger.error(str(e))
                failed_features.append({'feature_id': i, 'error': str(e), 'envelope': (minX, maxX, minY, maxY)})
                i += 1
                continue

            except MemoryError as e:
                error_msg = f"Memory error during GDAL Warp for feature {i}: {str(e)}"
                logger.error(error_msg)
                failed_features.append({'feature_id': i, 'error': error_msg, 'envelope': (minX, maxX, minY, maxY)})
                # Force garbage collection
                import gc
                gc.collect()
                i += 1
                continue

            except Exception as e:
                error_msg = f"Unexpected exception during GDAL Warp for feature {i}: {str(e)}"
                logger.error(error_msg)
                logger.error(f"Polygon envelope: {minX}, {maxX}, {minY}, {maxY}")
                logger.error(f"Exception type: {type(e).__name__}")
                logger.error(f"Traceback: {traceback.format_exc()}")
                # Log polygon geometry for debugging (truncated)
                try:
                    pPolygonWKT = pPolygon.ExportToWkt()
                    logger.debug(f"Polygon WKT (first 200 chars): {pPolygonWKT[:200]}...")
                except:
                    logger.debug("Could not log polygon WKT")
                failed_features.append({'feature_id': i, 'error': error_msg, 'envelope': (minX, maxX, minY, maxY)})
                i += 1
                continue

            # Optimize data type conversion and missing value handling
            if aData_clip is not None:
                # Use numpy's faster operations
                aData_clip = aData_clip.astype(np.int32, copy=False)  # Avoid unnecessary copy
                np.place(aData_clip, aData_clip == dMissing_value, -9999)  # Faster than boolean indexing

                if iFlag_save_clipped_raster_in == 1 :
                    if geometry_type == ogr.wkbPolygon : #we cannot save the multipolygon directly
                        sFilename_raster_out = os.path.join(sFolder_raster_out_in, f'{sRasterName_no_extension}_clip_{sClip}{sExtension}')
                        # Delete existing file if it exists (GDAL Create() doesn't overwrite)
                        if os.path.exists(sFilename_raster_out):
                            try:
                                pDriver.Delete(sFilename_raster_out)
                            except:
                                os.remove(sFilename_raster_out)  # Fallback if GDAL delete fails
                        iNewWidth = aData_clip.shape[1]
                        iNewHeigh = aData_clip.shape[0]
                        pDataset_clip = pDriver.Create(sFilename_raster_out, iNewWidth, iNewHeigh, 1, eType , options= options)
                        pDataset_clip.SetGeoTransform( newGeoTransform )
                        pDataset_clip.SetProjection( pSpatialRef_target_wkt)
                        #set the no data value
                        pBand = pDataset_clip.GetRasterBand(1)
                        pBand.SetNoDataValue(-9999)
                        pBand.WriteArray(aData_clip)
                        pBand.FlushCache()  # Corrected method name to FlushCache()
                        pDataset_clip.FlushCache()
                        pDataset_clip = None

                aCoords_gcs = get_geometry_coordinates(pPolygon)
                dArea = calculate_polygon_area(aCoords_gcs[:,0], aCoords_gcs[:,1])
                #create a polygon feature to save the output
                pFeature_out.SetGeometry(pPolygon.Clone())
                pFeature_out.SetField('id', i)
                pFeature_out.SetField('area', dArea)

                if iFlag_stat_in == 1:
                    # Check if aData_clip is already 1D valid data (from IDL merge) or 2D raster
                    if aData_clip.ndim == 1:
                        # Already 1D array of valid data from IDL merge
                        valid_data = aData_clip[aData_clip != -9999]  # Remove any remaining -9999 values
                    else:
                        # Standard 2D raster - extract valid data
                        valid_mask = aData_clip != -9999
                        valid_data = aData_clip[valid_mask]

                    if valid_data.size > 0:
                        # Calculate all statistics in one pass for efficiency
                        valid_data_float = valid_data.astype(np.float64, copy=False)
                        pFeature_out.SetField('mean', float(np.mean(valid_data_float)))
                        pFeature_out.SetField('min', float(np.min(valid_data_float)))
                        pFeature_out.SetField('max', float(np.max(valid_data_float)))
                        pFeature_out.SetField('std', float(np.std(valid_data_float)))

                pLayer_out.CreateFeature(pFeature_out)

                # Explicit cleanup - handle both single polygon dataset and multipolygon flag
                if isinstance(pDataset_clip_warped, bool):
                    # Multipolygon case - cleanup was already handled in helper method
                    pass
                else:
                    # Single polygon case - cleanup the dataset
                    pDataset_clip_warped = None

                # Progress reporting
                if i % 100 == 0:
                    elapsed = time.time() - start_time
                    rate = i / elapsed
                    eta = (len(aFeatures) - i) / rate if rate > 0 else 0
                    print(f"Processed {i}/{len(aFeatures)} features ({rate:.2f} features/sec, ETA: {eta:.0f}s)")

            i += 1


        pDataset_out.FlushCache()  # Flush once after all features are added
        pDataset_out = None        # Close the dataset
        pDataset_mesh = None

        # Report processing summary
        total_time = time.time() - start_time
        logger.info(f"Processing completed in {total_time:.2f} seconds")
        logger.info(f"Successfully processed: {successful_features} features")
        logger.info(f"Failed features: {len(failed_features)}")

        if failed_features:
            logger.warning("Failed features summary:")
            for failed in failed_features[:10]:  # Show first 10 failures
                logger.warning(f"  Feature {failed['feature_id']}: {failed['error']}")
            if len(failed_features) > 10:
                logger.warning(f"  ... and {len(failed_features) - 10} more failures")

        # Save failure report to file
        if failed_features:
            # Generate failure report filename by replacing extension with '_failures.log'
            base_name = os.path.splitext(sFilename_vector_out)[0]
            failure_report_file = f"{base_name}_failures.log"
            try:
                with open(failure_report_file, 'w') as f:
                    f.write(f"Processing failure report - {time.ctime()}\n")
                    f.write(f"Total features processed: {len(aFeatures)}\n")
                    f.write(f"Successful: {successful_features}\n")
                    f.write(f"Failed: {len(failed_features)}\n\n")
                    for failed in failed_features:
                        f.write(f"Feature {failed['feature_id']}: {failed['error']}\n")
                        f.write(f"  Envelope: {failed['envelope']}\n\n")
                logger.info(f"Failure report saved to: {failure_report_file}")
            except Exception as e:
                logger.error(f"Could not save failure report: {e}")

        return

    def _process_single_polygon(self, polygon, aFilename_source_raster, gdal_warp_options_base, feature_id):
        """
        Process a single polygon with GDAL Warp
        Returns: (dataset, data_array, geotransform) or (None, None, None) on failure
        """
        try:
            # Use WKT for faster performance
            pPolygonWKT = polygon.ExportToWkt()
            warp_options = gdal_warp_options_base.copy()
            warp_options['cutlineWKT'] = pPolygonWKT
            pWrapOption = gdal.WarpOptions(**warp_options)

            # Capture GDAL errors during Warp operation
            gdal.PushErrorHandler('CPLQuietErrorHandler')

            # Set a timeout for the operation (using alarm signal on Unix)
            def timeout_handler(signum, frame):
                raise TimeoutError(f"GDAL Warp operation timed out after 30 seconds for feature {feature_id}")

            if hasattr(signal, 'SIGALRM'):  # Unix systems only
                signal.signal(signal.SIGALRM, timeout_handler)
                signal.alarm(30)  # 30 second timeout

            logger.debug(f"Starting GDAL Warp for single polygon feature {feature_id}")
            pDataset_clip_warped = gdal.Warp('', aFilename_source_raster, options=pWrapOption)

            if hasattr(signal, 'SIGALRM'):
                signal.alarm(0)  # Cancel the alarm

            gdal.PopErrorHandler()

            # Check if the operation succeeded
            if pDataset_clip_warped is None:
                gdal_error = gdal.GetLastErrorMsg()
                logger.error(f"GDAL Warp failed for single polygon feature {feature_id}: {gdal_error}")
                return None, None, None

            newGeoTransform = pDataset_clip_warped.GetGeoTransform()
            aData_clip = pDataset_clip_warped.ReadAsArray()

            # Check if data was successfully read
            if aData_clip is None:
                logger.error(f"Failed to read array data for single polygon feature {feature_id}")
                pDataset_clip_warped = None
                return None, None, None

            # Check for reasonable data dimensions
            if aData_clip.size == 0:
                logger.warning(f"Empty data array for single polygon feature {feature_id}")
                pDataset_clip_warped = None
                return None, None, None

            return pDataset_clip_warped, aData_clip, newGeoTransform

        except TimeoutError as e:
            logger.error(str(e))
            return None, None, None
        except Exception as e:
            logger.error(f"Unexpected exception during single polygon processing for feature {feature_id}: {str(e)}")
            return None, None, None

    def _process_multipolygon_idl(self, multipolygon, aFilename_source_raster, gdal_warp_options_base, feature_id, dMissing_value):
        """
        Process a multipolygon (IDL crossing) by handling each part separately and merging results
        Returns: (merged_data_array, merged_geotransform) or (None, None) on failure
        """
        try:
            nGeometries = multipolygon.GetGeometryCount()
            logger.info(f"Processing {nGeometries} polygon parts for IDL-crossing feature {feature_id}")

            merged_datasets = []
            merged_data_arrays = []
            merged_transforms = []

            # Process each polygon part separately
            for iPart in range(nGeometries):
                polygon_part = multipolygon.GetGeometryRef(iPart)

                # Process this polygon part
                pDataset_part, aData_part, transform_part = self._process_single_polygon(
                    polygon_part, aFilename_source_raster, gdal_warp_options_base,
                    f"{feature_id}_part{iPart}")

                if pDataset_part is None:
                    logger.warning(f"Failed to process polygon part {iPart} of feature {feature_id}")
                    continue

                merged_datasets.append(pDataset_part)
                merged_data_arrays.append(aData_part)
                merged_transforms.append(transform_part)

            if not merged_data_arrays:
                logger.error(f"No polygon parts could be processed for IDL feature {feature_id}")
                return None, None

            # Merge the data arrays and transforms
            merged_data, merged_transform = self._merge_raster_parts(merged_data_arrays, merged_transforms, feature_id, dMissing_value)

            # Clean up datasets
            for dataset in merged_datasets:
                dataset = None

            return merged_data, merged_transform

        except Exception as e:
            logger.error(f"Error processing multipolygon IDL feature {feature_id}: {str(e)}")
            return None, None

    def _merge_raster_parts(self, data_arrays, transforms, feature_id, dMissing_value):
        """
        Merge multiple raster arrays from IDL-split polygons
        Returns only valid data as a 1D array for statistics calculation
        Returns: (merged_1D_array, dummy_transform)
        """
        try:
            if len(data_arrays) == 1:
                # Single array: extract valid data as 1D
                data_array = data_arrays[0]
                valid_mask = data_array != dMissing_value
                valid_data = data_array[valid_mask]
                return valid_data.flatten() if valid_data.size > 0 else np.array([dMissing_value]), transforms[0]

            # Multiple arrays: concatenate all valid data values
            all_valid_data = []

            for data_array in data_arrays:
                # Extract all valid (non-nodata) values
                valid_mask = data_array != dMissing_value
                valid_data = data_array[valid_mask]

                if valid_data.size > 0:
                    all_valid_data.append(valid_data.flatten())

            if not all_valid_data:
                logger.warning(f"No valid data found in any part for IDL feature {feature_id}")
                return np.array([dMissing_value]), transforms[0]

            # Concatenate all valid data into a single 1D array
            merged_data = np.concatenate(all_valid_data)

            logger.info(f"Successfully merged {len(data_arrays)} raster parts for feature {feature_id}: "
                       f"{merged_data.size} valid pixels")

            # Return 1D array of valid data only (no transform needed for statistics)
            return merged_data, transforms[0]

        except Exception as e:
            logger.error(f"Error merging raster parts for feature {feature_id}: {str(e)}")
            return None, None