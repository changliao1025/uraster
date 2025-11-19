import os, sys, platform
# Get the directory of the current script
sPath_current = os.path.dirname(os.path.abspath(__file__))
sPath_library = os.path.dirname(os.path.dirname(sPath_current))
sys.path.append(sPath_library)

# Construct the relative path to the data folder
sFolder_data = os.path.join(sPath_current, '..', '..', 'data', 'example_9')

# Print or use the data folder path
print(f"Data folder path: {sFolder_data}")

# Convert absolute paths to relative paths
sFilename_source_mesh = os.path.join(sFolder_data,  'input','tin.geojson') # use the L10-100 test mesh
sFilename_hydrosheds_dem = os.path.join(sFolder_data,  'input', 'hyd_glo_dem_15s.tif')

sFilename_target_mesh = os.path.join(sFolder_data, 'output', 'tin_uraster.geojson')
sFilename_mesh_png = os.path.join(sFolder_data,'output', 'mesh.jpg')
sFilename_variable_png = os.path.join(sFolder_data,'output', 'uraster.png')
sFilename_variable_animation = os.path.join(sFolder_data,'output', 'global_uraster.mp4')

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
    #use delaware area for visualization
    dLongitude_focus_in = -76.2033964
    dLatitude_focus_in = 39.739215
    pRaster.visualize_source_mesh(sFilename_out=sFilename_mesh_png,
                                   dLongitude_focus_in=dLongitude_focus_in,
                                     dLatitude_focus_in=dLatitude_focus_in,
                                     dZoom_factor = 0.7)

    pRaster.run_remap()
    #pRaster.report_outputs()
    sColormap = 'terrain'

    pRaster.visualize_target_mesh(
        sFilename_out=sFilename_variable_png,
        dLongitude_focus_in=dLongitude_focus_in,
        dZoom_factor = 7,
        dLatitude_focus_in=dLatitude_focus_in,
        sColormap=sColormap)

    #pRaster.visualize_target_mesh(
    #    sFilename_out=sFilename_variable_animation,
    #    sColormap=sColormap,
    #    dLongitude_focus_in=dLongitude_focus_in,
    #    dLatitude_focus_in=dLatitude_focus_in,
    #    iFlag_create_animation=True,
    #    iAnimation_frames=360,       # 1° longitude per frame
    #    sAnimation_format='mp4')

    #pRaster.cleanup()

    print('done')

if __name__ == '__main__':
    main()