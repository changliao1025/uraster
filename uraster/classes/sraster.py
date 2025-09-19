import os
import os
# Import GDAL spatial reference
from osgeo import osr
from osgeo import gdal
# Define a class named 'sraster'


class sraster:
	def __init__(self):
		# File path
		self.sFilepath = None

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
		self.spatial_ref = osr.SpatialReference()

		# Affine transform (geotransform)
		self.pTransform = None
		self.aExtent = None
		self.aExtent_wgs84 = None

		self.sFilename_mesh = None

		# NoData value
		self.dNoData = None
		return

	def read_metadata(self, sFilename):
		"""
		Read raster metadata from the given filename using GDAL.
		"""

		#check if file exists
		if not os.path.isfile(sFilename):
			raise FileNotFoundError(f"File does not exist: {sFilename}")

		pDataset = gdal.Open(sFilename)
		if pDataset is None:
			raise FileNotFoundError(f"Unable to open file: {sFilename}")

		self.sFilepath = sFilename
		self.iWidth = pDataset.RasterXSize
		self.iHeight = pDataset.RasterYSize
		self.iBandCount = pDataset.RasterCount
		self.sDtype = gdal.GetDataTypeName(pDataset.GetRasterBand(1).DataType) if self.iBandCount > 0 else None
		self.pTransform = pDataset.GetGeoTransform()
		self.dNoData = pDataset.GetRasterBand(1).GetNoDataValue() if self.iBandCount > 0 else None
		self.sCrs = pDataset.GetProjection()
		self.spatial_ref.ImportFromWkt(self.sCrs)
		pDataset = None

		#obtain the spatial extent
		self.aExtent = (self.pTransform[0], self.pTransform[3], self.pTransform[1], self.pTransform[4])

		#if the spatial reference is not in WGS84, transform the extent to WGS84
		if not self.spatial_ref.IsSame(osr.SRS_WKT_WGS84):
			transform = osr.CoordinateTransformation(self.spatial_ref, osr.SpatialReference().ImportFromEPSG(4326))
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
        """
        if not os.path.isfile(sFilename_mesh):
            raise FileNotFoundError(f"Mesh file does not exist: {sFilename_mesh}")

        self.sFilename_mesh = sFilename_mesh
        # Placeholder for actual mesh creation logic


