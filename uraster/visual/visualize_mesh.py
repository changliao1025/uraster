import os, sys
import numpy as np
from osgeo import ogr
import geovista as gv
#import geovista.theme
import geovista.report as gvreport
from pyearth.system.define_global_variables import *
from pyearth.gis.gdal.read.vector.get_supported_formats import get_format_from_extension
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.extract_unique_vertices_and_connectivity import extract_unique_vertices_and_connectivity

def visualize_mesh(sFilename_mesh, sFilename_out = None, sVariable_in = None, sUnit_in = None,
                    dLongitude_focus_in=0.0, dLatitude_focus_in=0.0, iFlag_show_gpu_info= None):

    if iFlag_show_gpu_info is not None:
        #check GPU availability
        print(gvreport.Report())

    if dLongitude_focus_in is not None:
        dLongitude_focus = dLongitude_focus_in
    else:
        dLongitude_focus = 0.0
    if dLatitude_focus_in is not None:
        dLatitude_focus = dLatitude_focus_in
    else:
        dLatitude_focus = 0.0

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

    if sVariable_in is None:
        #default to the first field
        sVariable = pLayerDefn.GetFieldDefn(0).GetName()
    else:
        sVariable = sVariable_in

    #loop the features to get lon/lat
    lons_list = []
    lats_list = []
    data_list = []
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
                lons_list.append(aCoord[:,0])
                lats_list.append(aCoord[:,1])
                data_list.append(pFeature.GetField(sVariable))


        pFeature = pLayer.GetNextFeature()

    # Pad to maximum size approach
    max_vertices = max(len(coord) for coord in lons_list)
    lons_padded = []
    lats_padded = []

    for lon_coords, lat_coords in zip(lons_list, lats_list):
        # Pad with NaN for missing vertices
        lon_padded = np.pad(lon_coords, (0, max_vertices - len(lon_coords)),
                           mode='constant', constant_values=np.nan)
        lat_padded = np.pad(lat_coords, (0, max_vertices - len(lat_coords)),
                           mode='constant', constant_values=np.nan)
        lons_padded.append(lon_padded)
        lats_padded.append(lat_padded)

    # Convert to numpy arrays
    lons = np.array(lons_padded)
    lats = np.array(lats_padded)

    # Get unique 1D coordinates for each mesh cell (centroids)
    cell_lons_1d = []
    cell_lats_1d = []

    for i in range(len(lons_list)):
        # Calculate centroid of each cell (ignoring NaN values)
        valid_lons = lons[i][~np.isnan(lons[i])]
        valid_lats = lats[i][~np.isnan(lats[i])]

        centroid_lon = np.mean(valid_lons)
        centroid_lat = np.mean(valid_lats)

        cell_lons_1d.append(centroid_lon)
        cell_lats_1d.append(centroid_lat)

    cell_lons_1d = np.array(cell_lons_1d)
    cell_lats_1d = np.array(cell_lats_1d)

    print(f"Number of mesh cells: {len(cell_lons_1d)}")
    print(f"Cell longitude range: {cell_lons_1d.min():.3f} to {cell_lons_1d.max():.3f}")
    print(f"Cell latitude range: {cell_lats_1d.min():.3f} to {cell_lats_1d.max():.3f}")

    # Extract unique vertices and connectivity using the independent function
    xv, yv, connectivity, vertex_to_index = extract_unique_vertices_and_connectivity(
        lons_list, lats_list
    )

    print(f"Number of unique vertices: {len(xv)}")
    print(f"Vertex longitude range: {xv.min():.3f} to {xv.max():.3f}")
    print(f"Vertex latitude range: {yv.min():.3f} to {yv.max():.3f}")
    print(f"Number of cells in connectivity: {len(connectivity)}")
    print(f"Sample connectivity for first cell: {connectivity[0] if len(connectivity) > 0 else 'None'}")

    ni, nj = lons.shape
    data = np.array(data_list)

    name = sVariable.capitalize()
    sUnit = sUnit_in if sUnit_in is not None else "unknown"

    # Create the mesh from the unstructured data using unique vertices and connectivity
    # Use masked connectivity for variable-sized polygons
    connectivity_masked = np.ma.masked_where(connectivity == -1, connectivity)
    crs = "EPSG:4326"
    mesh = gv.Transform.from_unstructured(xv, yv, connectivity=connectivity_masked, crs = crs)

    # Add the data values to the mesh
    print(f"Type of mesh.cell_data: {type(mesh.cell_data)}")
    print(f"MRO of mesh.cell_data: {type(mesh.cell_data).__mro__}")
    print(f"Data dtype before assignment: {data.dtype}")
    mesh.cell_data[name] = data
    print(f"Data dtype after assignment: {mesh.cell_data[name].dtype if name in mesh.cell_data.keys() else 'Not found'}")
    # Plot the mesh.
    plotter = gv.GeoPlotter()
    sargs = {"title": f"{name} / {sUnit}",
           "shadow": True,    "title_font_size": 10,    "label_font_size": 10,    "fmt": "%.1f",
    }
    plotter.add_mesh(mesh, scalars=name, scalar_bar_args=sargs)
    focal_point = gv.geodesic.to_cartesian([dLongitude_focus], [dLatitude_focus])[0]
    camera_position = gv.geodesic.to_cartesian([dLongitude_focus], [dLatitude_focus], radius=gv.common.RADIUS *6)[0]
    plotter.camera.focal_point = focal_point
    plotter.camera.position = camera_position
    plotter.camera.zoom(1.4)
    plotter.add_coastlines()
    plotter.add_axes()
    plotter.add_graticule(show_labels=True)
    if sFilename_out is not None:
        #save the figure
        plotter.screenshot(sFilename_out)
        print(f'Saved screenshot to: {sFilename_out}')
    else:
        plotter.show()
