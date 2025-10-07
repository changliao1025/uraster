import os, sys, platform



sPath  = 'C:\\workspace\\python\\uraster\\uraster'
sys.path.append(os.path.dirname(sPath)) #add the parent path to the system path

from uraster.visual.visualize_mesh import visualize_mesh
sPlatform_os = platform.system()
if sPlatform_os == 'win32':
    sFilename_target_mesh = 'C:\\scratch\\04model\\pyhexwatershed\\global\\pyflowline20250926004\\mpas.json' #use the L10-100 test mesh
else:
    sFilename_target_mesh = 'Z:\\04model\\pyhexwatershed\\global\\pyflowline20250926004\\mpas.geojson' #use the L10-100 test mesh

visualize_mesh(sFilename_target_mesh)