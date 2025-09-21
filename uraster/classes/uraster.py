# Define a class named 'uraster'
import os
import numpy as np
from osgeo import gdal, ogr, osr
from .sraster import sraster

class uraster:

    def __init__(self):
        self.iFlag_global = None
        self.dResolution_raster = None
        self.dResolution_uraster = None

        self.sFilename_source_mesh = None
        self.sFilename_target_mesh = None


    def check_raster_files(self, aFilename_source_raster_in):
        """
        Check if the input raster files exist
        :param aFilename_source_raster_in: a list of raster files
        :return: True if all files exist, False otherwise
        """
        for sFilename_raster_in in aFilename_source_raster_in:
            if os.path.exists(sFilename_raster_in):
                pass
            else:
                print('The raster file does not exist:', sFilename_raster_in)
                return False
        #uset the sraster class the check the raster
        for sFilename_raster_in in aFilename_source_raster_in:
            pRaster = sraster(sFilename_raster_in)
            pRaster.read_metadata()
            



        return True



    def remap_raster_to_uraster(self, aFilename_source_raster_in,    sFilename_target_mesh_in ,
                  iFlag_save_clipped_raster_in = None,
                  iFlag_save_in_vector_in = None,
                  iFlag_stat_in = None,
                  iFlag_return_data_in = None,
                  sFilename_vector_out= None,
                  sFolder_raster_out= None,
                  sFormat_in='GTiff'):
        """
        Clip a raster by a shapefile
        :param aFilename_source_raster_in: a list of raster files
        :param    sFilename_target_mesh_in : input polygon filename
        :param sFilename_raster_out: output raster filename
        :param sFormat: output format
        :return: None
        """

        #check input files
        for sFilename_raster_in in aFilename_source_raster_in:
            if os.path.exists(sFilename_raster_in):
                pass
            else:
                print('The raster file does not exist!')
                return

        if os.path.exists(   sFilename_target_mesh_in ):
            pass
        else:
            print('The shapefile does not exist!')
            return

        pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')

        #check the input raster data format and decide gdal driver
        if sFormat_in is not None:
            sDriverName = sFormat_in
        else:
            sDriverName = 'GTiff'

        if iFlag_save_clipped_raster_in is not None:
            iFlag_save_clipped_raster = iFlag_save_clipped_raster_in
        else:
            iFlag_save_clipped_raster = 0

        if iFlag_save_in_vector_in is not None:
            iFlag_save_in_vector = iFlag_save_in_vector_in
        else:
            iFlag_save_in_vector = 0

        if iFlag_save_in_vector == 1:
            if os.path.exists(sFilename_vector_out):
                #remove the file using the gdal driver
                pDriver_shapefile.DeleteDataSource(sFilename_vector_out)

        if iFlag_save_clipped_raster == 1:
            if sFolder_raster_out is not None:
                if os.path.exists(sFolder_raster_out):
                    pass
                else:
                    #create the folder
                    os.makedirs(sFolder_raster_out)

            else:
                #save the output to the same folder as the input
                sFolder_raster_out = os.path.dirname(sFilename_raster_in)
        else:
            print('The output raster will not be saved unless there is only one polygon!')

        pDriver = gdal.GetDriverByName(sDriverName)

        sFilename_shapefile_cut = "/vsimem/tmp_polygon.shp"
        #get the raster file extension
        sExtension = os.path.splitext(sFilename_raster_in)[1]
        sName = os.path.basename(sFilename_raster_in)
        sRasterName_no_extension = os.path.splitext(sName)[0]

        pDataset_data = gdal.Open(sFilename_raster_in, gdal.GA_ReadOnly)

        dummy = gdal_read_geotiff_file(sFilename_raster_in)

        eType = dummy['dataType']
        dPixelWidth = dummy['pixelWidth']
        pPixelHeight = dummy['pixelHeight']
        dOriginX = dummy['originX']
        dOriginY = dummy['originY']
        nrow = dummy['nrow']
        ncolumn = dummy['ncolumn']
        dMissing_value= dummy['missingValue']
        pProjection = dummy['projection']
        pSpatialRef_target= osr.SpatialReference()
        pSpatialRef_target.ImportFromWkt(pProjection)
        wkt1 = pProjection
        # Check if the raster's spatial reference is projected
        is_raster_projected = pSpatialRef_target.IsProjected()
        is_raster_geographic = pSpatialRef_target.IsGeographic()

        dX_left=dOriginX
        dX_right = dOriginX + ncolumn * dPixelWidth
        dY_top = dOriginY
        dY_bot = dOriginY + nrow * pPixelHeight

        #get the spatial reference of the shapefile
        pDataset_subset = ogr.Open(   sFilename_target_mesh_in )
        pLayer = pDataset_subset.GetLayer(0)
        # Count the number of features (polygons)
        nFeature = pLayer.GetFeatureCount()
        # Get the spatial reference of the layer
        pSpatialRef_source = pLayer.GetSpatialRef()
        wkt2 = pSpatialRef_source.ExportToWkt()

        #check whether the polygon has only one or more features
        if nFeature > 1:
            pass
        else:
            print('The polygon file has only one polygons!')
            return

        #check the polygon spatial reference, reproject if necessary
        if(wkt1 != wkt2):
            pDataset_clip = None
            #in this case, we can reproject the shapefile to the same spatial reference as the raster
            #get the folder that contains the shapefile
            sFolder = os.path.dirname(sFilename_vector_out)
            #get the name of the shapefile
            sName = os.path.basename(   sFilename_target_mesh_in )
            #get the name of the shapefile without extension
            sName_no_extension = os.path.splitext(sName)[0]
            #create a new shapefile
            sFilename_clip_out = sFolder + '/' + sName_no_extension + '_transformed.shp'
            reproject_vector(   sFilename_target_mesh_in , sFilename_clip_out, pProjection)
            #use the new shapefile to clip the raster
            sFilename_clip = sFilename_clip_out
        else:
            sFilename_clip =    sFilename_target_mesh_in
            #read the first polygon

        #get the envelope of the polygon

        #now loop through all the polygons in the vector file
        pDataset_subset = None
        pLayer_subset = None
        pFeature_subset = None
        pDataset_subset = ogr.Open(sFilename_clip)
        pLayer_subset = pDataset_subset.GetLayer(0)
        aData_subset = list()
        options = ['COMPRESS=DEFLATE', 'PREDICTOR=2']
        sResampleAlg='near'
        sDriverName = 'MEM'
        pDriver_memory = gdal.GetDriverByName(sDriverName)
        if iFlag_save_clipped_raster == 1:
            print('The clipped raster will be saved to the folder:', sFolder_raster_out)
            pass

        if iFlag_save_in_vector == 1:
            #create a polygon feature to save the output
            pDataset2 = pDriver_shapefile.CreateDataSource(sFilename_vector_out)
            pLayerOut2 = pDataset2.CreateLayer('cell', pSpatialRef_target, ogr.wkbPolygon)
            pLayerDefn2 = pLayerOut2.GetLayerDefn()
            pFeatureOut2 = ogr.Feature(pLayerDefn2)
            #add id, area and mean, min, max, std of the raster
            pLayerOut2.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
            #define a field
            pField = ogr.FieldDefn('area', ogr.OFTReal)
            pField.SetWidth(32)
            pField.SetPrecision(2)
            pLayerOut2.CreateField(pField)
            if iFlag_stat_in == 1:
                pLayerOut2.CreateField(ogr.FieldDefn('mean', ogr.OFTReal))
                pLayerOut2.CreateField(ogr.FieldDefn('min', ogr.OFTReal))
                pLayerOut2.CreateField(ogr.FieldDefn('max', ogr.OFTReal))
                pLayerOut2.CreateField(ogr.FieldDefn('std', ogr.OFTReal))
            else:
                #treat the raster as categorical
                pLayerOut2.CreateField(ogr.FieldDefn('ntype', ogr.OFTInteger))
                pLayerOut2.CreateField(ogr.FieldDefn('type', ogr.OFTInteger))
                pLayerOut2.CreateField(ogr.FieldDefn('percent', ogr.OFTReal))
                pass
            pass

            #all the clipped rasters will be saved to the output folder
            #for i in range(nFeature):

        i = 1
        pLayer_subset.ResetReading()
        pFeature_subset = pLayer_subset.GetNextFeature()
        while pFeature_subset is not None:
            sClip = "{:06d}".format(i)
            pPolygon = pFeature_subset.GetGeometryRef()
            minX, maxX, minY, maxY = pPolygon.GetEnvelope()
            if minX > dX_right or maxX < dX_left  or minY > dY_top or maxY < dY_bot:
                pass
            else:
                if not pPolygon or pPolygon.IsEmpty() or not pPolygon.IsValid():
                    pFeature_subset = pLayer_subset.GetNextFeature()
                    continue
                #create the temporary shapefile
                #pDataset3 = pDriver_shapefile.CreateDataSource(sFilename_shapefile_cut)
                #pLayerOut3 = pDataset3.CreateLayer('cell', pSpatialRef_target, ogr.wkbPolygon)
                #pLayerDefn3 = pLayerOut3.GetLayerDefn()
                #pFeatureOut3 = ogr.Feature(pLayerDefn3)
                #pFeatureOut3.SetGeometry(pPolygon)
                #pLayerOut3.CreateFeature(pFeatureOut3)
                #pDataset3.FlushCache()
                #pWrapOption = gdal.WarpOptions( cropToCutline=True,cutlineDSName = sFilename_shapefile_cut , \
                #        width=iNewWidth, \
                #            height=iNewHeigh, \
                #                dstSRS=pProjection , format = sDriverName )
                pPolygonWKT = pPolygon.ExportToWkt()
                pWrapOption = gdal.WarpOptions(cropToCutline=True,
                                                #cutlineDSName = sFilename_shapefile_cut ,#could be true if vector file is provided
                                                cutlineWKT=pPolygonWKT,
                                    xRes=dPixelWidth,
                                   yRes=abs(pPixelHeight),
                                        dstSRS=pSpatialRef_target , format = 'MEM',
                                        resampleAlg=sResampleAlg )
                pDataset_clip_warped = gdal.Warp('', pDataset_data,
                                                  options=pWrapOption)#gdal.Warp(sFilename_out, pDataset_in, options=pWrapOption)
                newGeoTransform = pDataset_clip_warped.GetGeoTransform()
                aData_clip = pDataset_clip_warped.ReadAsArray().astype(int)
                aData_clip[aData_clip == dMissing_value] = -9999

                if iFlag_save_clipped_raster == 1:
                    sFilename_raster_out  = sFolder_raster_out + '/' + sRasterName_no_extension + '_clip_' + sClip + sExtension
                    iNewWidth = aData_clip.shape[1]
                    iNewHeigh = aData_clip.shape[0]
                    pDataset_clip = pDriver.Create(sFilename_raster_out, iNewWidth, iNewHeigh, 1, eType , options= options)
                    pDataset_clip.SetGeoTransform( newGeoTransform )
                    pDataset_clip.SetProjection( pProjection)
                    #set the no data value
                    pBand = pDataset_clip.GetRasterBand(1)
                    pBand.SetNoDataValue(-9999)
                    pBand.WriteArray(aData_clip)
                    pBand.FlushCache()  # Corrected method name to FlushCache()
                    pDataset_clip.FlushCache()
                    pDataset_clip = None

                if is_raster_projected:
                    dArea= pPolygon.GetArea()
                else:
                    aCoords_gcs = get_geometry_coordinates(pPolygon)
                    dArea = calculate_polygon_area(aCoords_gcs[:,0], aCoords_gcs[:,1])
                if iFlag_save_in_vector == 1:
                    #create a polygon feature to save the output
                    pFeatureOut2.SetGeometry(pPolygon.Clone())
                    pFeatureOut2.SetField('id', i)
                    pFeatureOut2.SetField('area', dArea)
                    if iFlag_stat_in == 1:
                        dummy_data = aData_clip[np.where(aData_clip != -9999)].astype(float)
                        if len(dummy_data) == 0:
                            pFeature_subset = pLayer_subset.GetNextFeature()
                            continue
                        else:
                            pFeatureOut2.SetField('mean', np.mean(dummy_data))
                            pFeatureOut2.SetField('min', np.min(dummy_data))
                            pFeatureOut2.SetField('max', np.max(dummy_data))
                            pFeatureOut2.SetField('std', np.std(dummy_data))
                    else:
                        #treat the raster as categorical
                        aData_clip = aData_clip.flatten()
                        dummy_index = np.where(aData_clip != -9999)
                        if len(dummy_index) == 0:
                            #no data
                            pFeatureOut2.SetField('ntype', 0)
                            pFeatureOut2.SetField('type', -9999)
                            pFeatureOut2.SetField('percent', 0)
                            pFeature_subset = pLayer_subset.GetNextFeature()
                            continue
                        else:
                            dummy_data = aData_clip[dummy_index]
                            if len(dummy_data) == 0:
                                pFeatureOut2.SetField('ntype', 0)
                                pFeatureOut2.SetField('type', -9999)
                                pFeatureOut2.SetField('percent', 0)
                                pFeature_subset = pLayer_subset.GetNextFeature()
                                continue
                            else:
                                #find the dominant value
                                aUnique, aCount = np.unique(dummy_data, return_counts=True)
                                iIndex = np.argmax(aCount)
                                pFeatureOut2.SetField('ntype', len(aUnique))
                                pFeatureOut2.SetField('type', int(aUnique[iIndex]))
                                #get percentage of the dominant value
                                iTotal = len(dummy_data)
                                iCount = aCount[iIndex]
                                p = float(iCount)/float(iTotal)
                                pFeatureOut2.SetField('percent',p)

                        pass
                    pLayerOut2.CreateFeature(pFeatureOut2)

                if iFlag_return_data_in is not None:
                    aData_subset.append(aData_clip)
                pDataset_clip_warped = None

            i = i + 1
            pFeature_subset = pLayer_subset.GetNextFeature()

        # After the loop finishes processing all features
        if iFlag_save_in_vector == 1 and pDataset2 is not None:
            pDataset2.FlushCache()  # Flush once after all features are added
            pDataset2 = None        # Close the dataset

        # ... (rest of your function, like closing other datasets) ...
        pDataset_data = None
        pDataset_subset = None

        return