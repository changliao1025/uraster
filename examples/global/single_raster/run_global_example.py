import os, sys, platform
sPlatform_os = platform.system()

if sPlatform_os == 'Windows':
    sPath  = 'C:\\workspace\\python\\uraster\\uraster'
    sys.path.append(os.path.dirname(sPath))
    sFilename_source_mesh = 'C:\\scratch\\04model\\pyhexwatershed\\global\\pyflowline20250927006\\mpas.geojson' #use the L10-100 test mesh
    sFilename_hydrosheds_dem = 'C:\\scratch\\00raw\\hydrology\\hydrosheds\\hydrosheds\\hyd_glo_dem_15s.tif'
    sFilename_target_mesh = 'C:\\scratch\\04model\\pyhexwatershed\\global\\pyflowline20250927006\\mpas_uraster.geojson'
    sFilename_mesh_png = 'C:\\workspace\\python\\uraster\\data\\output\\mpas\\mesh.jpg'
    sFilename_variable_png = 'C:\\workspace\\python\\uraster\\data\\output\\mpas\\uraster.png'
    sFilename_variable_animation = 'C:\\workspace\\python\\uraster\\data\\output\\mpas\\global_uraster.mp4'
else:
    #macOS
    if sPlatform_os == 'Darwin':
        sFilename_source_mesh = '/Users/liao313/scratch/04model/pyhexwatershed/global/pyflowline20250927006//mpas.geojson' #use the L10-100 test mesh
        sFilename_hydrosheds_dem = '/Users/liao313/scratch/00raw/hydrology/hydrosheds/hydrosheds/hyd_glo_dem_15s.tif'
        sFilename_target_mesh = '/Users/liao313/scratch/04model/pyhexwatershed/global/pyflowline20250927006/mpas_uraster.geojson'
    else:
        #linux
        sFilename_source_mesh = '/compyfs/liao313/04model/pyhexwatershed/global/pyflowline20250927006/mpas.geojson' #use the L10-100 test mesh
        sFilename_hydrosheds_dem = '/compyfs/liao313/00raw/hydrology/hydrosheds/hydrosheds/hyd_glo_dem_15s.tif'
        sFilename_target_mesh = '/compyfs/liao313/04model/pyhexwatershed/global/pyflowline20250927006/mpas_uraster.geojson'

from uraster.classes.uraster import uraster

def main():
    aConfig=dict()
    aConfig['sFilename_source_mesh']= sFilename_source_mesh #use the L10-100 test mesh
    aFilename_source_raster = []

    aFilename_source_raster.append(sFilename_hydrosheds_dem) #dem from hydros
    aConfig['aFilename_source_raster']= aFilename_source_raster
    aConfig['sFilename_target_mesh']= sFilename_target_mesh
    pRaster = uraster(aConfig)

    pRaster.setup()
    pRaster.report_inputs()
    dLongitude_focus_in = -112.033964
    dLatitude_focus_in = 43.491977
    pRaster.visualize_source_mesh(sFilename_out=sFilename_mesh_png, dLongitude_focus_in=dLongitude_focus_in, dLatitude_focus_in=dLatitude_focus_in)

    pRaster.run_remap()
    #pRaster.report_outputs()
    sColormap = 'terrain'

    pRaster.visualize_target_mesh(
        sFilename_out=sFilename_variable_png,
        dLongitude_focus_in=dLongitude_focus_in,
        dLatitude_focus_in=dLatitude_focus_in,
        sColormap=sColormap)

    pRaster.visualize_target_mesh(
        sFilename_out=sFilename_variable_animation,
        sColormap=sColormap,
        dLongitude_focus_in=dLongitude_focus_in,
        dLatitude_focus_in=dLatitude_focus_in,
        iFlag_create_animation=True,
        iAnimation_frames=360,       # 1° longitude per frame
        sAnimation_format='mp4')

    #pRaster.cleanup()

    print('done')

if __name__ == '__main__':
    main()