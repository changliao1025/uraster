import os, sys, platform
sPath_current = os.path.dirname(os.path.abspath(__file__))
sPath_library = os.path.dirname(os.path.dirname(sPath_current))
sys.path.append(sPath_library)

# Construct the relative path to the data folder
sFolder_data = os.path.join(sPath_current, '..', '..', 'data', 'example_1')
sFolder_data = os.path.realpath(sFolder_data)
# Print or use the data folder path
print(f"Data folder path: {sFolder_data}")


# Convert absolute paths to relative paths
sFilename_source_mesh = os.path.join(sFolder_data, 'input',  'rhealpix_global_res3.geojson') # use the L10-100 test mesh
sFilename_raster = os.path.join(sFolder_data, 'input', 'EDGAR_CH4_emission_global_2015.tiff')

sFilename_target_mesh = os.path.join(sFolder_data, 'output',  'uraster.geojson')
sFilename_mesh_png = os.path.join(sFolder_data, 'output',  'mesh.jpg')
sFilename_raster_png = os.path.join(sFolder_data, 'output',  'raster.png')
sFilename_variable_png = os.path.join(sFolder_data, 'output',  'uraster.png')
sFilename_variable_animation = os.path.join(sFolder_data, 'output', 'uraster.gif')


from uraster.classes.uraster import uraster

def main():
    aConfig = dict()
    aConfig['sFilename_source_mesh'] = sFilename_source_mesh  # use the L10-100 test mesh
    aFilename_source_raster = []

    aFilename_source_raster.append(sFilename_raster)  # dem from hydros
    aConfig['aFilename_source_raster'] = aFilename_source_raster
    aConfig['sFilename_target_mesh'] = sFilename_target_mesh
    pRaster = uraster(aConfig)

    pRaster.setup()

    pRaster.report_inputs()
    # visualize source mesh at the Idaho Falls area
    dLongitude_focus_in = (pRaster.aExtent_rasters[0] + pRaster.aExtent_rasters[1]) / 2
    dLatitude_focus_in = (pRaster.aExtent_rasters[2] + pRaster.aExtent_rasters[3]) / 2
    pRaster.visualize_source_mesh(sFilename_out=sFilename_mesh_png, dLongitude_focus_in=dLongitude_focus_in, dLatitude_focus_in=dLatitude_focus_in)
    #pRaster.visualize_raster(sFilename_out=sFilename_raster_png)

    pRaster.run_remap()
    pRaster.report_outputs() #not implemented yet
    sColormap = 'terrain'

    #Optional visualization and animation (disabled by default in this script)
    pRaster.visualize_target_mesh(
        sFilename_out=sFilename_variable_png,
        sColormap=sColormap, dLongitude_focus_in=dLongitude_focus_in, dLatitude_focus_in=dLatitude_focus_in)

    pRaster.visualize_target_mesh(
        sFilename_out=sFilename_variable_animation,
        sColormap=sColormap,
        dLongitude_focus_in=dLongitude_focus_in, dLatitude_focus_in=dLatitude_focus_in,
        iFlag_create_animation=True,
        iAnimation_frames=360,       # 1° longitude per frame
        sAnimation_format='mp4')

    pRaster.cleanup()

    print('done')


if __name__ == '__main__':
    main()