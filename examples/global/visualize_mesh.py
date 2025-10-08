import os, sys, platform

sPath  = 'C:\\workspace\\python\\uraster\\uraster'
sys.path.append(os.path.dirname(sPath))

from uraster.visual.visualize_mesh import visualize_mesh
sPlatform_os = platform.system()
if sPlatform_os == 'Windows':
    sFilename_target_mesh = 'C:\\scratch\\04model\\pyhexwatershed\\global\\pyflowline20250926004\\mpas.geojson' #use the L10-100 test mesh
else:
    sFilename_target_mesh = 'Z:\\04model\\pyhexwatershed\\global\\pyflowline20250926004\\mpas.geojson' #use the L10-100 test mesh

visualize_mesh(sFilename_target_mesh, dLongitude_focus_in=-179.0, dLatitude_focus_in=67.0, iFlag_show_gpu_info=1)