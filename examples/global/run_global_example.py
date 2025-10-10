import os, sys, platform
sPlatform_os = platform.system()

if sPlatform_os == 'Windows':
    sPath  = 'C:\\workspace\\python\\uraster\\uraster'
    sys.path.append(os.path.dirname(sPath))
    sFilename_target_mesh = 'C:\\scratch\\04model\\pyhexwatershed\\global\\pyflowline20250926004\\mpas.geojson' #use the L10-100 test mesh
else:
    sFilename_target_mesh = 'Z:\\04model\\pyhexwatershed\\global\\pyflowline20250926004\\mpas.geojson' #use the L10-100 test mesh

from uraster.classes.uraster import uraster

aConfig=dict()

aConfig['sFilename_target_mesh']= sFilename_target_mesh #use the L10-100 test mesh
aFilename_source_raster = []
sFilename_hydrosheds_dem = '/compyfs/liao313/00raw/hydrology/hydrosheds/hydrosheds/hyd_glo_dem_15s.tif'
aFilename_source_raster.append(sFilename_hydrosheds_dem) #dem from hydros
pRaster = uraster(aConfig)

pRaster.setup()
pRaster.report_inputs()
pRaster.visualize_target_mesh()
sFilename_out_vector = '/compyfs/liao313/04model/pyhexwatershed/global/pyflowline20250926003/mpas_uraster.geojson'
pRaster.run_remap(sFilename_out_vector)
pRaster.report_outputs()
pRaster.visualize_outputs()