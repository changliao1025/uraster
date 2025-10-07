import os
import glob
from uraster.classes.uraster import uraster


#create a uraster object using the dict

aConfig=dict()
sFilename_target_mesh = '/compyfs/liao313/04model/pyhexwatershed/global/pyflowline20250926003/mpas.geojson' #use the L10-100 test mesh
aConfig['sFilename_target_mesh']= sFilename_target_mesh #use the L10-100 test mesh

aFilename_source_raster = list()
#find all the tif under the gebco folder
sFolder_gebco = '/compyfs/liao313/00raw/dem/global/gebco/normal/'
sPattern_gebco = os.path.join(sFolder_gebco, '*.tif')
for sFilename in glob.glob(sPattern_gebco):
    aFilename_source_raster.append(sFilename) #dem from gebco 2020

aConfig['aFilename_source_raster']=aFilename_source_raster  #resolution of the input raster in degrees


sFilename_hydrosheds_dem = '/compyfs/liao313/00raw/hydrology/hydrosheds/hydrosheds/hyd_glo_dem_15s.tif'
aFilename_source_raster.clear()
aFilename_source_raster.append(sFilename_hydrosheds_dem) #dem from hydros
pRaster = uraster(aConfig)

#print the uraster object attributes

pRaster.check_raster_files()
pRaster.print_raster_info()
#specify the output file, which is a vector file similar to the mesh file
sFilename_out_vector = '/compyfs/liao313/04model/pyhexwatershed/global/pyflowline20250926003/mpas_uraster.shp'

#sFolder_raster_out = '/compyfs/liao313/04model/pyhexwatershed/global/dem'
pRaster.remap_raster_to_uraster(sFilename_out_vector)