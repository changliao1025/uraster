import os
import os
# Import GDAL spatial reference
from osgeo import osr
from osgeo import gdal
# Define a class named 'sraster'
from pyearth.toolbox.management.raster.reproject import reproject_raster_gdalwarp


class sraster:
    def __init__(self, sFilename=None):
        # File path
        self.sFilename = sFilename

        # Raster dimensions
        self.iWidth = None
        self.iHeight = None

        # Number of bands
        self.iBandCount = None

        # Data type (e.g., uint8, float32)
        self.sDtype = None

        # Coordinate Reference System (CRS)
        self.sCrs = None

        # GDAL Spatial Reference object
        self.pSpatialRef = osr.SpatialReference()
        self.pSpatialRef_wkt = None

        # Affine transform (geotransform)
        self.pTransform = None
        self.aExtent = None
        self.aExtent_wgs84 = None

        self.sFilename_mesh = None

        # NoData value
        self.dNoData = None
        return

    def read_metadata(self):
        """
        Read raster metadata from the given filename using GDAL.
        """

        #check if file exists
        sFilename = self.sFilename
        if not os.path.isfile(sFilename):
            raise FileNotFoundError(f"File does not exist: {sFilename}")

        pDataset = gdal.Open(sFilename)
        if pDataset is None:
            raise FileNotFoundError(f"Unable to open file: {sFilename}")

        self.sFilename = sFilename
        self.iWidth = pDataset.RasterXSize
        self.iHeight = pDataset.RasterYSize
        self.iBandCount = pDataset.RasterCount
        self.sDtype = gdal.GetDataTypeName(pDataset.GetRasterBand(1).DataType) if self.iBandCount > 0 else None
        self.pTransform = pDataset.GetGeoTransform()
        self.dNoData = pDataset.GetRasterBand(1).GetNoDataValue() if self.iBandCount > 0 else None
        self.sCrs = pDataset.GetProjection()
        self.pSpatialRef.ImportFromWkt(self.sCrs)
        self.pSpatialRef_wkt = self.pSpatialRef.ExportToWkt() if self.sCrs else None
        pDataset = None

        #obtain the spatial extent
        self.aExtent = (self.pTransform[0], self.pTransform[3], self.pTransform[1], self.pTransform[4])

        #if the spatial reference is not in WGS84, transform the extent to WGS84
        if not self.pSpatialRef.IsSame(osr.SRS_WKT_WGS84):
            transform = osr.CoordinateTransformation(self.pSpatialRef, osr.SpatialReference().ImportFromEPSG(4326))
            (minX, maxY, _) = transform.TransformPoint(self.aExtent[0], self.aExtent[1])
            (maxX, minY, _) = transform.TransformPoint(self.aExtent[0] + self.aExtent[2], self.aExtent[1] + self.aExtent[3])
            self.aExtent_wgs84 = (minX, minY, maxX, maxY)
        else:
            self.aExtent_wgs84 = self.aExtent

        return

    def create_raster_mesh(self):
        """
        Create a raster mesh from the given mesh filename.
        """

        #check raster file sFilename exists
        if not os.path.isfile(self.sFilename):
            raise FileNotFoundError(f"Raster file does not exist: {self.sFilename}")

        #check mesh file exists, if yes, delete it
        if self.sFilename_mesh and os.path.isfile(self.sFilename_mesh):
            os.remove(self.sFilename_mesh)

        # Placeholder for actual mesh creation logic
        #we need to use pyflowline approach to create the mesh because?

    def convert_to_wgs84(self):
        """
        Convert the raster to WGS84 coordinate system.
        """
        # Define the target spatial reference (WGS84)
        pSpatialRef_wgs84 = osr.SpatialReference()
        pSpatialRef_wgs84.ImportFromEPSG(4326)
        pSpatialRef_wgs84_wkt = pSpatialRef_wgs84.ExportToWkt()

        #define a new filename for the converted raster
        sFilename_raster_wgs84 = self.sFilename.replace('.tif', '_wgs84.tif')
        #delete the file if it already exists
        if os.path.isfile(sFilename_raster_wgs84):
            os.remove(sFilename_raster_wgs84)
        #use/copy a function from the pyearth package to do the conversion
        reproject_raster_gdalwarp(self.sFilename, sFilename_raster_wgs84, pSpatialRef_wgs84_wkt,
                              xRes=None, yRes=None,
                               sResampleAlg = 'near',
                                 iFlag_force_resolution_in = 0)

        self.sFilename_wgs84 = sFilename_raster_wgs84

        return sraster(sFilename=sFilename_raster_wgs84)






