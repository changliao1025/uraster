import os, sys, platform

from pyearth.toolbox.management.raster.resample import resample_raster

sPath_current = os.path.dirname(os.path.abspath(__file__))
sPath_library = os.path.dirname(os.path.dirname(sPath_current))
sys.path.append(sPath_library)
from uraster.classes.uraster import uraster
# Download input data using Pooch (downloads to system cache)
from uraster.utility import get_example_paths

print("Downloading example 6 input data...")
paths = get_example_paths(example_number=6)
sFolder_input = paths['input']
print(f"Input data cached at: {sFolder_input}")

# Set up output directory relative to current working directory
sFolder_output = os.path.join("data", "example_6", "output")
os.makedirs(sFolder_output, exist_ok=True)
print(f"Output directory: {sFolder_output}")

# Convert absolute paths to relative paths
sFilename_source_mesh = os.path.join(
    sFolder_input, "isea7h_bbox_res8.geojson"
)  # use the L10-100 test mesh

sFilename_target_mesh = os.path.join(sFolder_output, "uraster.geojson")
sFilename_mesh_png = os.path.join(sFolder_output, "mesh.jpg")
sFilename_raster_png = os.path.join(sFolder_output, "raster.png")
sFilename_variable_png = os.path.join(sFolder_output, "uraster.png")
sFilename_variable_animation = os.path.join(sFolder_output, "uraster.mp4")




def main():
    aConfig = dict()
    aConfig["sFilename_source_mesh"] = (
        sFilename_source_mesh  # use the L10-100 test mesh
    )
    aFilename_source_raster = []
    sFilename_raster0 = os.path.join(
        sFolder_input, "ArcticDEM_32_34_1_1_2m_v4.1_dem.tif"
    )
    sFilename_raster1 = os.path.join(
        sFolder_input, "ArcticDEM_32_34_1_2_2m_v4.1_dem.tif"
    )
    sFilename_raster2 = os.path.join(
        sFolder_input, "ArcticDEM_32_35_1_1_2m_v4.1_dem.tif"
    )
    sFilename_raster3 = os.path.join(sFolder_input, "N79W078_FABDEM_V1-0.tif")
    sFilename_raster4 = os.path.join(sFolder_input, "N79W079_FABDEM_V1-0.tif")
    sFilename_raster5 = os.path.join(sFolder_input, "N79W080_FABDEM_V1-0.tif")
    sFilename_raster6 = os.path.join(sFolder_input, "N79W081_FABDEM_V1-0.tif")
    sFilename_raster7 = os.path.join(sFolder_input, "N79W082_FABDEM_V1-0.tif")
    sFilename_raster8 = os.path.join(sFolder_input, "N79W083_FABDEM_V1-0.tif")
    aFilename_source_raster.append(sFilename_raster0)  #
    aFilename_source_raster.append(sFilename_raster1)  #
    aFilename_source_raster.append(sFilename_raster2)  #
    aFilename_source_raster.append(sFilename_raster3)  #
    aFilename_source_raster.append(sFilename_raster4)  #
    aFilename_source_raster.append(sFilename_raster5)  #
    aFilename_source_raster.append(sFilename_raster6)  #
    aFilename_source_raster.append(sFilename_raster7)  #
    aFilename_source_raster.append(sFilename_raster8)  #
    aConfig["aFilename_source_raster"] = aFilename_source_raster
    aConfig["sFilename_target_mesh"] = sFilename_target_mesh
    # use weighted average remap method
    pRaster = uraster(aConfig)
    pRaster.setup()
    pRaster.report_inputs()
    # visualize source mesh at the Wuhan City area
    dLongitude_focus_in = (pRaster.aExtent_rasters[0] + pRaster.aExtent_rasters[2]) / 2
    dLatitude_focus_in = (pRaster.aExtent_rasters[1] + pRaster.aExtent_rasters[3]) / 2
    pRaster.visualize_source_mesh(
        sFilename_out=sFilename_mesh_png,
        dLongitude_focus_in=dLongitude_focus_in,
        dLatitude_focus_in=dLatitude_focus_in,
    )
    # pRaster.visualize_raster(sFilename_out=sFilename_raster_png)

    pRaster.run_remap(iFlag_weighted_average_in=False)
    pRaster.report_outputs()  # not implemented yet
    sColormap = "terrain"
    # Optional visualization and animation (disabled by default in this script)
    pRaster.visualize_target_mesh(
        sFilename_out=sFilename_variable_png,
        sColormap=sColormap,
        dLongitude_focus_in=dLongitude_focus_in,
        dLatitude_focus_in=dLatitude_focus_in,
    )

    pRaster.visualize_target_mesh(
        sFilename_out=sFilename_variable_animation,
        sColormap=sColormap,
        dLongitude_focus_in=dLongitude_focus_in,
        dLatitude_focus_in=dLatitude_focus_in,
        iFlag_create_animation=True,
        iAnimation_frames=360,  # 1° longitude per frame
        sAnimation_format="mp4",
    )

    pRaster.cleanup()
    print("done")


if __name__ == "__main__":
    main()
