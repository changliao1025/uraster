import os, sys, platform

# Get the directory of the current script
sPath_current = os.path.dirname(os.path.abspath(__file__))
sPath_library = os.path.dirname(os.path.dirname(sPath_current))
sys.path.append(sPath_library)
from uraster.classes.uraster import uraster
# Download input data using Pooch (downloads to system cache)
from uraster.utility import get_example_paths

print("Downloading example 8 input data...")
paths = get_example_paths(example_number=8)
sFolder_input = paths['input']
print(f"Input data cached at: {sFolder_input}")

# Set up output directory relative to current working directory
sFolder_output = os.path.join("data", "example_8", "output")
os.makedirs(sFolder_output, exist_ok=True)
print(f"Output directory: {sFolder_output}")

# Convert absolute paths to relative paths
sFilename_source_mesh = os.path.join(
    sFolder_input, "tin.geojson"
)  # use the L10-100 test mesh
sFilename_hydrosheds_dem = os.path.join(sFolder_input, "hyd_glo_dem_15s.tif")

sFilename_target_mesh = os.path.join(sFolder_output, "uraster.geojson")
sFilename_mesh_png = os.path.join(sFolder_output, "mesh.jpg")
sFilename_variable_png = os.path.join(sFolder_output, "uraster.png")
sFilename_variable_animation = os.path.join(sFolder_output, "uraster.mp4")




def main():
    aConfig = dict()
    aConfig["sFilename_source_mesh"] = (
        sFilename_source_mesh  # use the L10-100 test mesh
    )
    aFilename_source_raster = []

    aFilename_source_raster.append(sFilename_hydrosheds_dem)  # dem from hydros
    aConfig["aFilename_source_raster"] = aFilename_source_raster
    aConfig["sFilename_target_mesh"] = sFilename_target_mesh
    pRaster = uraster(aConfig)
    pRaster.setup()
    pRaster.report_inputs()
    # use delaware area for visualization
    dLongitude_focus_in = -75.0
    dLatitude_focus_in = 39.0
    pRaster.visualize_source_mesh(
        sFilename_out=sFilename_mesh_png,
        dLongitude_focus_in=dLongitude_focus_in,
        dLatitude_focus_in=dLatitude_focus_in,
        dZoom_factor=0.7,
    )

    pRaster.run_remap()
    pRaster.report_outputs()
    sColormap = "terrain"

    pRaster.visualize_target_mesh(
        sFilename_out=sFilename_variable_png,
        dLongitude_focus_in=dLongitude_focus_in,
        dLatitude_focus_in=dLatitude_focus_in,
        sColormap=sColormap,
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
