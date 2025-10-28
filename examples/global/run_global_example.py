import os, sys, platform
sPlatform_os = platform.system()

if sPlatform_os == 'Windows':
    sPath  = 'C:\\workspace\\python\\uraster\\uraster'
    sys.path.append(os.path.dirname(sPath))
    sFilename_source_mesh = 'C:\\scratch\\04model\\pyhexwatershed\\global\\pyflowline20250926008\\mpas.geojson' #use the L10-100 test mesh
    sFilename_hydrosheds_dem = 'Z:\\00raw\\hydrology\\hydrosheds\\hydrosheds\\hyd_glo_dem_15s.tif'
    sFilename_out_vector = 'C:\\scratch\\04model\\pyhexwatershed\\global\\pyflowline20250927006\\mpas_uraster.geojson'
else:
    #macOS
    if sPlatform_os == 'Darwin':
        sFilename_source_mesh = '/Users/liao313/scratch/04model/pyhexwatershed/global/pyflowline20250926008//mpas.geojson' #use the L10-100 test mesh
        sFilename_hydrosheds_dem = '/Users/liao313/scratch/00raw/hydrology/hydrosheds/hydrosheds/hyd_glo_dem_15s.tif'
        sFilename_out_vector = '/Users/liao313/scratch/04model/pyhexwatershed/global/pyflowline20250927006/mpas_uraster.geojson'
    else:
        #linux
        sFilename_source_mesh = '/compyfs/liao313/04model/pyhexwatershed/global/pyflowline20250926008/mpas.geojson' #use the L10-100 test mesh
        sFilename_hydrosheds_dem = '/compyfs/liao313/00raw/hydrology/hydrosheds/hydrosheds/hyd_glo_dem_15s.tif'
        sFilename_out_vector = '/compyfs/liao313/04model/pyhexwatershed/global/pyflowline20250927006/mpas_uraster.geojson'

from uraster.classes.uraster import uraster
aConfig=dict()
aConfig['sFilename_source_mesh']= sFilename_source_mesh #use the L10-100 test mesh
aFilename_source_raster = []

aFilename_source_raster.append(sFilename_hydrosheds_dem) #dem from hydros
aConfig['aFilename_source_raster']= aFilename_source_raster
aConfig['sFilename_target_mesh']= sFilename_out_vector
pRaster = uraster(aConfig)

pRaster.setup()
#pRaster.report_inputs()
#pRaster.visualize_source_mesh()

#pRaster.run_remap(sFilename_out_vector)
#pRaster.report_outputs()
pRaster.visualize_target_mesh()

print('done')