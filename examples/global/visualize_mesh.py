import os
import glob
from uraster.visual.visualize_mesh import visualize_mesh


sFilename_target_mesh = '/compyfs/liao313/04model/pyhexwatershed/global/pyflowline20250926003/mpas.geojson' #use the L10-100 test mesh

visualize_mesh(sFilename_target_mesh)