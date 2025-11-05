import os, sys, platform
sPlatform_os = platform.system()


if sPlatform_os == 'Darwin':
    sFilename_source_mesh = '/Users/liao313/workspace/python/uraster/data/input/mesh/gdf_rhpx_res3.geojson' #use the L10-100 test mesh
    sFilename_hydrosheds_dem = '/Users/liao313/workspace/python/uraster/data/input/raster/edgar_MNM_2015.tiff'
    sFilename_target_vector = '/Users/liao313/workspace/python/uraster/data/output/rhpx/rhpx_res3_uraster.geojson'


from uraster.classes.uraster import uraster
aConfig=dict()
aConfig['sFilename_source_mesh']= sFilename_source_mesh #use the L10-100 test mesh
aFilename_source_raster = []

aFilename_source_raster.append(sFilename_hydrosheds_dem) #dem from hydros
aConfig['aFilename_source_raster']= aFilename_source_raster
aConfig['sFilename_target_mesh']= sFilename_target_vector
pRaster = uraster(aConfig)

pRaster.setup()
pRaster.report_inputs()
pRaster.visualize_source_mesh(sFilename_out='/Users/liao313/workspace/python/uraster/data/output/rhpx/mesh.jpg')
exit()
pRaster.run_remap()
pRaster.report_outputs()
sColormap = 'terrain'

pRaster.visualize_target_mesh(
    sFilename_out='/Users/liao313/workspace/python/uraster/data/output/rhpx/uraster.png',
    sColormap=sColormap)

pRaster.visualize_target_mesh(
    sFilename_out='/Users/liao313/workspace/python/uraster/data/output/rhpx/global_uraster.mp4',
    sColormap=sColormap,
    iFlag_create_animation=True,
    iAnimation_frames=360,       # 1° longitude per frame
    sAnimation_format='mp4')

pRaster.cleanup()

print('done')