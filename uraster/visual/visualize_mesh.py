import os, sys
import numpy as np
from osgeo import ogr
import geovista as gv
import geovista.theme
from pyearth.gis.gdal.read.vector.get_supported_formats import get_format_from_extension
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
def visualize_mesh(sFilename_mesh, sFilename_out = None, sVariable_in = None, sUnit_in = None):

    #call gdal api to get the lon/lat from the vector file

    #determine file extension
    sExt = os.path.splitext(sFilename_mesh)[1]
    sFormat = get_format_from_extension(sFilename_mesh)
    print(f'Auto-detected format: {sFormat}')

    # Get the appropriate driver
    pDriver_out = ogr.GetDriverByName(sFormat)
    if pDriver_out is None:
        print(f'{sFormat} driver not available.')
        return

    # Open the input data source
    pDataset = ogr.Open(sFilename_mesh, 0)  # 0 means read-only. 1 means writeable.
    if pDataset is None:
        print(f'Failed to open file: {sFilename_mesh}')
        return
    print(f'Opened file: {sFilename_mesh}')
    # Get the first layer
    pLayer = pDataset.GetLayer()
    if pLayer is None:
        print('Failed to get layer from the dataset.')
        return

    # Get the layer definition
    pLayerDefn = pLayer.GetLayerDefn()
    if pLayerDefn is None:
        print('Failed to get layer definition.')
        return
    iFieldCount = pLayerDefn.GetFieldCount()
    print(f'Number of fields in the layer: {iFieldCount}')

    #loop the features to get lon/lat
    lons = []
    lats = []
    #reset the reading to the start
    pLayer.ResetReading()
    pFeature = pLayer.GetNextFeature()
    while pFeature is not None:
        pGeometry = pFeature.GetGeometryRef()
        if pGeometry is not None:
            # Assuming the geometry is a polygon, get the centroid
            sGeometry_type = pGeometry.GetGeometryName()
            if sGeometry_type == 'POLYGON':
                #get all the coordinates of the polygon
                aCoord = get_geometry_coordinates(pGeometry)
                lons.append(aCoord[:,0])
                lats.append(aCoord[:,1])


        pFeature = pLayer.GetNextFeature()


    #reshape to remove the second dimension

    lons = np.array(lons)
    lats = np.array(lats)
    ni, nj, nv = lons.shape
    lons = lons.reshape(ni, nv)
    lats = lats.reshape(ni, nv)
    connectivity = lons.shape

    if sVariable_in is None:
        #default to the first field
        sVariable = pLayerDefn.GetFieldDefn(0).GetName()
    else:
        sVariable = sVariable_in

    name = sVariable.capitalize()
    sUnit = sUnit_in if sUnit_in is not None else "unknown"

    # Create the mesh from the sample data.
    t = 0
    mesh = gv.Transform.from_unstructured(lons, lats)
    # Plot the mesh.
    plotter = gv.GeoPlotter()
    sargs = {"title": f"{name} / {sUnit}",
           "shadow": True,    "title_font_size": 10,    "label_font_size": 10,    "fmt": "%.1f",
    }
    plotter.add_mesh(mesh, scalar_bar_args=sargs)
    plotter.view_xy()


    plotter.camera.zoom(1.4)
    plotter.add_coastlines()
    plotter.add_axes()
    if sFilename_out is not None:
        #save the figure
        plotter.screenshot(sFilename_out)
        print(f'Saved screenshot to: {sFilename_out}')
    else:
        plotter.show()
