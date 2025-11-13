import os, sys, platform
sPlatform_os = platform.system()

if sPlatform_os == 'Windows':
    sPath  = 'C:\\workspace\\python\\uraster\\uraster'
    sys.path.append(os.path.dirname(sPath))
    sFilename_source_mesh = 'C:\\workspace\\python\\uraster\\data\\input\\mesh\\gdf_rhpx_res3.geojson' #use the L10-100 test mesh
    sFilename_raster = 'C:\\workspace\\python\\uraster\\data\\input\\raster\\edgar_MNM_2015.tiff'
    sFilename_target_mesh = 'C:\\workspace\\python\\uraster\\data\\output\\rhpx\\rhpx_res3_uraster.geojson'
    sFilename_mesh_png = 'C:\\workspace\\python\\uraster\\data\\output\\rhpx\\mesh.jpg'
    sFilename_variable_png = 'C:\\workspace\\python\\uraster\\data\\output\\rhpx\\uraster.png'
    sFilename_variable_animation = 'C:\\workspace\\python\\uraster\\data\\output\\rhpx\\global_uraster.mp4'
else:
    if sPlatform_os == 'Darwin':
        sFilename_source_mesh = '/Users/liao313/workspace/python/uraster/data/input/mesh/gdf_rhpx_res3.geojson' #use the L10-100 test mesh
        sFilename_raster = '/Users/liao313/workspace/python/uraster/data/input/raster/edgar_MNM_2015.tiff'
        sFilename_target_mesh = '/Users/liao313/workspace/python/uraster/data/output/rhpx/rhpx_res3_uraster.geojson'
        sFilename_mesh_png = '/Users/liao313/workspace/python/uraster/data/output/rhpx/mesh.jpg'
        sFilename_raster_png = '/Users/liao313/workspace/python/uraster/data/output/rhpx/raster.png'
        sFilename_variable_png = '/Users/liao313/workspace/python/uraster/data/output/rhpx/uraster.png'
        sFilename_variable_animation = '/Users/liao313/workspace/python/uraster/data/output/rhpx/global_uraster.mp4'


from uraster.classes.uraster import uraster

def main():
    aConfig = dict()
    aConfig['sFilename_source_mesh'] = sFilename_source_mesh  # use the L10-100 test mesh
    aFilename_source_raster = []

    aFilename_source_raster.append(sFilename_raster)  # dem from hydros
    aConfig['aFilename_source_raster'] = aFilename_source_raster
    aConfig['sFilename_target_mesh'] = sFilename_target_mesh
    pRaster = uraster(aConfig)

    pRaster.setup()

    pRaster.report_inputs()
    # visualize source mesh at the Idaho Falls area
    dLongitude_focus_in = -112.033964
    dLatitude_focus_in = 43.491977
    pRaster.visualize_source_mesh(sFilename_out=sFilename_mesh_png, dLongitude_focus_in=dLongitude_focus_in, dLatitude_focus_in=dLatitude_focus_in)
    #pRaster.visualize_raster(sFilename_out=sFilename_raster_png)
    #exit()
    pRaster.run_remap()
    #pRaster.report_outputs() #not implemented yet
    sColormap = 'terrain'

    #Optional visualization and animation (disabled by default in this script)
    pRaster.visualize_target_mesh(
        sFilename_out=sFilename_variable_png,
        sColormap=sColormap, dLongitude_focus_in=dLongitude_focus_in, dLatitude_focus_in=dLatitude_focus_in)

    pRaster.visualize_target_mesh(
        sFilename_out=sFilename_variable_animation,
        sColormap=sColormap,
        dLongitude_focus_in=dLongitude_focus_in, dLatitude_focus_in=dLatitude_focus_in,
        iFlag_create_animation=True,
        iAnimation_frames=360,       # 1° longitude per frame
        sAnimation_format='mp4')

    pRaster.cleanup()

    print('done')


if __name__ == '__main__':
    main()