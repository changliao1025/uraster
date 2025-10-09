from uraster.classes.uraster import uraster

aConfig=dict()
sFilename_target_mesh = '/compyfs/liao313/04model/pyhexwatershed/global/pyflowline20250926003/mpas.geojson' #use the L10-100 test mesh
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