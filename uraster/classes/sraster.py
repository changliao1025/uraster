import os
import numpy as np
from osgeo import osr
from osgeo import gdal
gdal.UseExceptions()
from pyearth.toolbox.management.raster.reproject import reproject_raster
from pyearth.toolbox.mesh.square.create_square_mesh import create_square_mesh
from pyearth.toolbox.mesh.latlon.create_latlon_mesh import create_latlon_mesh

class sraster:

    def __init__(self, sFilename_in=None):
        # File path

        self.sFilename = sFilename_in
        # Raster dimensions
        self.ncolumn = None
        self.nrow = None
        # Number of bands
        self.iBandCount = None
        # Data type (e.g., uint8, float32)
        self.sDtype = None
        # Coordinate Reference System (CRS)
        self.sCrs = None
        # GDAL Spatial Reference object - initialize to prevent AttributeError
        self.pSpatialRef = osr.SpatialReference()
        self.pSpatialRef_wkt = None

        # Affine transform (geotransform)
        self.pTransform = None
        self.aExtent = None
        self.aExtent_wgs84 = None
        self.sFilename_mesh = None
        #check the file exists
        if sFilename_in is not None:
            if os.path.isfile(sFilename_in):
                #setup the mesh filename using robust path handling
                base, ext = os.path.splitext(sFilename_in)
                self.sFilename_mesh = f"{base}_mesh.geoparquet"
            else:
                raise FileNotFoundError(f"File does not exist: {sFilename_in}")
        # NoData value
        self.dNoData = None

        # Resolution attributes
        self.dResolution_x = None
        self.dResolution_y = None
        self.nrow = None
        self.ncolumn = None

        # WGS84 extent bounds
        self.dLongitude_left = None
        self.dLongitude_right = None
        self.dLatitude_bottom = None
        self.dLatitude_top = None
        return

    def read_metadata(self):
        """
        Read raster metadata from the given filename using GDAL.
        """
        pSpatialRef_wgs84 = osr.SpatialReference()
        pSpatialRef_wgs84.ImportFromEPSG(4326)
        wkt_wgs84 = pSpatialRef_wgs84.ExportToWkt()
        pSpatialRef_wgs84 = None  # Clean up
        #check if file exists
        sFilename = self.sFilename
        if not os.path.isfile(sFilename):
            raise FileNotFoundError(f"File does not exist: {sFilename}")

        pDataset = None
        try:
            pDataset = gdal.Open(sFilename)
            if pDataset is None:
                raise RuntimeError(f"Unable to open file: {sFilename}")

            self.sCrs = pDataset.GetProjection()
            if self.sCrs:
                self.pSpatialRef.ImportFromWkt(self.sCrs)
                self.pSpatialRef_wkt = self.pSpatialRef.ExportToWkt()
            else:
                self.pSpatialRef_wkt = None

            self.ncolumn = pDataset.RasterXSize
            self.nrow = pDataset.RasterYSize
            self.iBandCount = pDataset.RasterCount
            self.eType = pDataset.GetRasterBand(1).DataType if self.iBandCount > 0 else None
            self.sDtype = gdal.GetDataTypeName(self.eType) if self.iBandCount > 0 else None
            self.pTransform = pDataset.GetGeoTransform()
            self.dResolution_x = self.pTransform[1]
            self.dResolution_y = -self.pTransform[5]
            self.dNoData = pDataset.GetRasterBand(1).GetNoDataValue() if self.iBandCount > 0 else None

            # Calculate the actual spatial extent (minX, minY, maxX, maxY)
            minX = self.pTransform[0]
            maxY = self.pTransform[3]
            maxX = minX + (self.ncolumn * self.pTransform[1])
            minY = maxY + (self.nrow * self.pTransform[5])  # pTransform[5] is negative
            self.aExtent = (minX, minY, maxX, maxY)

            # If the spatial reference is not in WGS84, transform the extent to WGS84
            if self.pSpatialRef_wkt != wkt_wgs84:
                pSpatialRef_wgs84_target = None
                transform = None
                try:
                    pSpatialRef_wgs84_target = osr.SpatialReference()
                    pSpatialRef_wgs84_target.ImportFromEPSG(4326)
                    transform = osr.CoordinateTransformation(self.pSpatialRef, pSpatialRef_wgs84_target)

                    # Transform all 4 corners to handle rotated rasters
                    (x1, y1, _) = transform.TransformPoint(minX, minY)  # Lower-left
                    (x2, y2, _) = transform.TransformPoint(maxX, minY)  # Lower-right
                    (x3, y3, _) = transform.TransformPoint(maxX, maxY)  # Upper-right
                    (x4, y4, _) = transform.TransformPoint(minX, maxY)  # Upper-left

                    # Get the bounding box of all transformed corners
                    all_x = [x1, x2, x3, x4]
                    all_y = [y1, y2, y3, y4]
                    self.aExtent_wgs84 = (min(all_x), min(all_y), max(all_x), max(all_y))

                except Exception as e:
                    print(f"Warning: Failed to transform extent to WGS84: {e}")
                    self.aExtent_wgs84 = self.aExtent
                finally:
                    # Clean up spatial reference and transformation objects
                    if pSpatialRef_wgs84_target is not None:
                        pSpatialRef_wgs84_target = None
                    if transform is not None:
                        transform = None
            else:
                self.aExtent_wgs84 = self.aExtent

            self.dLongitude_left = self.aExtent_wgs84[0]
            self.dLongitude_right = self.aExtent_wgs84[2]
            self.dLatitude_bottom = self.aExtent_wgs84[1]
            self.dLatitude_top = self.aExtent_wgs84[3]

        finally:
            # Ensure dataset is properly closed
            if pDataset is not None:
                pDataset = None

        return

    def print_info(self):
        """
        Print raster metadata information.
        """
        print(f"Filename: {self.sFilename}")
        print(f"Width: {self.ncolumn}")
        print(f"Height: {self.nrow}")
        print(f"Band Count: {self.iBandCount}")
        print(f"Data Type: {self.sDtype}")
        print(f"NoData Value: {self.dNoData}")
        print(f"CRS: {self.sCrs}")
        print(f"Spatial Reference WKT: {self.pSpatialRef_wkt}")
        print(f"Affine Transform: {self.pTransform}")
        print(f"Extent: {self.aExtent}")
        print(f"WGS84 Extent: {self.aExtent_wgs84}")
        return

    def create_raster_mesh(self):
        """
        Create a raster mesh from the raster file.
        Creates either a square mesh (for projected CRS) or a lat/lon mesh (for WGS84).
        """
        pSpatialRef_wgs84 = osr.SpatialReference()
        pSpatialRef_wgs84.ImportFromEPSG(4326)
        wkt_wgs84 = pSpatialRef_wgs84.ExportToWkt()
        pSpatialRef_wgs84 = None  # Clean up

        #check raster file sFilename exists
        if not os.path.isfile(self.sFilename):
            raise FileNotFoundError(f"Raster file does not exist: {self.sFilename}")

        # Ensure metadata has been read
        if self.ncolumn is None or self.nrow is None:
            raise RuntimeError("Metadata not loaded. Call read_metadata() first.")

        #check mesh file exists, if yes, delete it
        if self.sFilename_mesh and os.path.isfile(self.sFilename_mesh):
            os.remove(self.sFilename_mesh)

        # Create mesh based on coordinate system
        if self.pSpatialRef_wkt != wkt_wgs84:
            # Projected coordinate system - use square mesh
            dX_left_in = self.aExtent[0]
            dY_bot_in = self.aExtent[1]
            dResolution_meter_in = self.dResolution_x  # Fixed: removed trailing comma
            ncolumn_in = self.ncolumn
            nrow_in = self.nrow
            sFilename_output_in = self.sFilename_mesh
            pProjection_reference_in = self.pSpatialRef_wkt

            create_square_mesh(dX_left_in, dY_bot_in,
                        dResolution_meter_in,
                        ncolumn_in, nrow_in,
                        sFilename_output_in,
                        pProjection_reference_in)
        else:
            # WGS84 geographic coordinate system - use lat/lon mesh
            dLongitude_left_in = self.dLongitude_left
            if dLongitude_left_in < -180.0:
                dLongitude_left_in =-180
            dLatitude_bot_in = self.dLatitude_bottom
            if dLatitude_bot_in < -90.0:
                dLatitude_bot_in = -90.0
            dResolution_degree_in = self.dResolution_x
            ncolumn_in = self.ncolumn
            nrow_in = self.nrow
            sFilename_output_in = self.sFilename_mesh

            create_latlon_mesh(dLongitude_left_in,
                       dLatitude_bot_in,
                       dResolution_degree_in,
                       ncolumn_in, nrow_in,
                       sFilename_output_in)

        return self.sFilename_mesh

    def convert_to_wgs84(self):
        """
        Convert the raster to WGS84 coordinate system.
        Returns a new sraster instance for the reprojected file.
        """
        pSpatialRef_wgs84 = osr.SpatialReference()
        pSpatialRef_wgs84.ImportFromEPSG(4326)
        wkt_wgs84 = pSpatialRef_wgs84.ExportToWkt()
        pSpatialRef_wgs84 = None  # Clean up
        # Check if already in WGS84
        if self.pSpatialRef_wkt == wkt_wgs84:
            print("Raster is already in WGS84 coordinate system.")
            return self


        # Define a new filename for the converted raster using robust path handling
        base, ext = os.path.splitext(self.sFilename)
        sFilename_raster_wgs84 = f"{base}_wgs84{ext}"

        # Delete the file if it already exists
        if os.path.isfile(sFilename_raster_wgs84):
            os.remove(sFilename_raster_wgs84)

        # Use/copy a function from the pyearth package to do the conversion
        reproject_raster(self.sFilename, sFilename_raster_wgs84, wkt_wgs84,
                              xRes = None, yRes = None,
                               sResampleAlg = 'near',
                                 iFlag_force_resolution_in = 0)

        # Create and return a new sraster instance for the reprojected file
        pRaster_wgs84 = sraster(sFilename_in=sFilename_raster_wgs84)
        pRaster_wgs84.read_metadata()
        return pRaster_wgs84

    def get_unique_values(self, band_index=1, dMissing_value=None, iFlag_verbose_in=0):
        """
        Get unique values from the specified band of the raster.
        """
        # Check if file exists
        if not os.path.isfile(self.sFilename):
            raise FileNotFoundError(f"File does not exist: {self.sFilename}")

        pDataset = None
        unique_values = set()
        try:
            pDataset = gdal.Open(self.sFilename)
            if pDataset is None:
                raise RuntimeError(f"Unable to open file: {self.sFilename}")

            if band_index < 1 or band_index > pDataset.RasterCount:
                raise ValueError(f"Invalid band index: {band_index}")

            pBand = pDataset.GetRasterBand(band_index)
            array = pBand.ReadAsArray()

            # Get unique values using numpy
            unique_array = np.unique(array)

            # Filter out missing value if specified
            for value in unique_array:
                if dMissing_value is not None and value == dMissing_value:
                    continue
                unique_values.add(value)

        finally:
            # Ensure dataset is properly closed
            if pDataset is not None:
                pDataset = None

        return list(unique_values)




