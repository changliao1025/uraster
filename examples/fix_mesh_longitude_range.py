import logging

# Set up logging
logger = logging.getLogger(__name__)
from osgeo import ogr, osr
crs = "EPSG:4326"
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename
from pyearth.gis.geometry.international_date_line_utility import split_international_date_line_polygon_coordinates, check_cross_international_date_line_polygon
def fix_idl_crossing(sFilename_mesh, sFilename_mesh_fixed):
    pDriver = get_vector_driver_from_filename(sFilename_mesh)
    #get the source dataset
    pDataset = pDriver.Open(sFilename_mesh, 0)  # 0 means
    if pDataset is None:
        logger.error(f"Could not open file: {sFilename_mesh}")
        return None
    pLayer = pDataset.GetLayer()
    pLayerDefn = pLayer.GetLayerDefn()
    if pLayer is None:
        logger.error(f"Could not get layer from file: {sFilename_mesh}")
        pDataset = None
        return None
    nFeatures = pLayer.GetFeatureCount()
    sSpatial_ref = pLayer.GetSpatialRef()

    #create the target dataset
    pDriver_out = get_vector_driver_from_filename(sFilename_mesh_fixed)
    pDataset_out = pDriver_out.CreateDataSource(sFilename_mesh_fixed)
    if pDataset_out is None:
        logger.error(f"Could not create file: {sFilename_mesh_fixed}")
        return None
    pLayer_out = pDataset_out.CreateLayer('layer', sSpatial_ref, ogr.wkbPolygon)
    if pLayer_out is None:
        logger.error(f"Could not create layer in file: {sFilename_mesh_fixed}")
        pDataset_out = None
        return None

    for pFeature in pLayer:
        geometry = pFeature.GetGeometryRef()
        geometry_type = geometry.GetGeometryName()
        if geometry_type == 'POLYGON':
            aCoord = get_geometry_coordinates(geometry)
            bCross_idl = check_cross_international_date_line_polygon(aCoord)
            if bCross_idl:
                logger.info(f"Feature ID {pFeature.GetFID()} crosses the International Date Line. Splitting...")
                [eastern_polygon, western_polygon] = split_international_date_line_polygon_coordinates(aCoord)
                #create a multiplepolygon geometry

                pGeometery_multi = ogr.Geometry(ogr.wkbMultiPolygon)
                #create eastern polygon feature
                pPolygon_eastern = ogr.Geometry(ogr.wkbPolygon)
                pLinearRing_eastern = ogr.Geometry(ogr.wkbLinearRing)
                for coord in eastern_polygon:
                    pLinearRing_eastern.AddPoint(coord[0], coord[1])
                pLinearRing_eastern.CloseRings()
                pPolygon_eastern.AddGeometry(pLinearRing_eastern)
                pGeometery_multi.AddGeometry(pPolygon_eastern)
                #create western polygon feature
                pPolygon_western = ogr.Geometry(ogr.wkbPolygon)
                pLinearRing_western = ogr.Geometry(ogr.wkbLinearRing)
                for coord in western_polygon:
                    pLinearRing_western.AddPoint(coord[0], coord[1])
                pLinearRing_western.CloseRings()
                pPolygon_western.AddGeometry(pLinearRing_western)
                pGeometery_multi.AddGeometry(pPolygon_western)

                #save this feature
                pFeature_out = ogr.Feature(pLayer_out.GetLayerDefn())
                pFeature_out.SetGeometry(pGeometery_multi)
                #copy over all fields
                for iField in range(pLayerDefn.GetFieldCount()):
                    sField_name = pLayerDefn.GetFieldDefn(iField).GetName()
                    pFeature_out.SetField(sField_name, pFeature.GetField(sField_name))
                pLayer_out.CreateFeature(pFeature_out)

            else:
                #save this feature directly
                pLayer_out.CreateFeature(pFeature)
        else:
            logger.warning(f"Geometry type {geometry_type} not supported for IDL crossing check.")


    #flush thedata to disk
    pDataset_out.FlushCache()
    pDataset_out = None
    pDataset = None


    pass