import os, sys, platform

sPath_current = os.path.dirname(os.path.abspath(__file__))
sPath_library = os.path.dirname(os.path.dirname(sPath_current))
sys.path.append(sPath_library)
from uraster.classes.uraster import uraster
# Download input data using Pooch (downloads to system cache)
from uraster.utility import get_example_paths

print("Downloading example 2 input data...")
paths = get_example_paths(example_number=2)
sFolder_input = paths['input']
print(f"Input data cached at: {sFolder_input}")

# Set up output directory relative to current working directory
sFolder_output = os.path.join("data", "example_2", "output")
os.makedirs(sFolder_output, exist_ok=True)
print(f"Output directory: {sFolder_output}")

# Input file paths (from Pooch cache)
sFilename_source_mesh = os.path.join(
    sFolder_input, "isea3h_bbox_res15.geojson"
)  # use the L10-100 test mesh
sFilename_raster = os.path.join(sFolder_input, "WSF2015_v2_-80_42.tif")

# Output file paths (relative to working directory)
sFilename_target_mesh = os.path.join(sFolder_output, "uraster.geojson")
sFilename_mesh_png = os.path.join(sFolder_output, "mesh.jpg")
sFilename_raster_png = os.path.join(sFolder_output, "raster.png")
sFilename_variable_png = os.path.join(sFolder_output, "uraster.png")
sFilename_variable_animation = os.path.join(sFolder_output, "uraster.gif")





def main():
    aConfig = dict()
    aConfig["sFilename_source_mesh"] = (
        sFilename_source_mesh  # use the L10-100 test mesh
    )
    aFilename_source_raster = []

    aFilename_source_raster.append(sFilename_raster)  # dem from hydros
    aConfig["aFilename_source_raster"] = aFilename_source_raster
    aConfig["sFilename_target_mesh"] = sFilename_target_mesh
    aConfig["iFlag_discrete"] = 1
    pRaster = uraster(aConfig)

    pRaster.setup()

    pRaster.report_inputs()
    # visualize source mesh at the Idaho Falls area
    dLongitude_focus_in = (pRaster.aExtent_rasters[0] + pRaster.aExtent_rasters[2]) / 2
    dLatitude_focus_in = (pRaster.aExtent_rasters[1] + pRaster.aExtent_rasters[3]) / 2
    pRaster.visualize_source_mesh(
        sFilename_out=sFilename_mesh_png,
        dLongitude_focus_in=dLongitude_focus_in,
        dLatitude_focus_in=dLatitude_focus_in,
    )
    # pRaster.visualize_raster(sFilename_out=sFilename_raster_png)

    pRaster.run_remap(iFlag_discrete_in=True)
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
