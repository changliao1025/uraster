import os

from osgeo import osr
from pyearth.toolbox.management.raster.resample import resample_raster

from uraster import URaster, compile_remap_workflow, load_workflow_config


def main():
    # 1. Load Tier 1 configuration from workflow_config.yaml
    sFilename_yaml = os.path.join(os.path.dirname(__file__), "workflow_config.yaml")
    config = load_workflow_config(sFilename_yaml)

    # Example 10 downsamples the raw input raster for the demonstration
    sFilename_raster = config["rasters"]["source_rasters"][0]
    sFilename_source_raster_resample = sFilename_raster.replace(".tif", "_resample.tif")
    if not os.path.exists(sFilename_source_raster_resample):
        default_srs = osr.SpatialReference()
        default_srs.ImportFromEPSG(4326)
        pProjection_target = default_srs.ExportToWkt()
        resample_raster(
            sFilename_raster,
            sFilename_source_raster_resample,
            1.0,
            1.0,
            sResampleAlg="average",
            dMissing_value_source=0,
            pProjection_target_in=pProjection_target,
        )

    # Update workflow config with resampled raster path
    config["rasters"]["source_rasters"] = [sFilename_source_raster_resample]

    # 2. Compile Tier 1 specification into backend execution parameters
    aConfig = compile_remap_workflow(config)
    pRaster = URaster(aConfig)

    # 3. Setup and report input metadata
    pRaster.setup()
    pRaster.report_inputs()

    # 4. Visualize source mesh
    pRaster.visualize_source_mesh()

    # 5. Run remap (mass-conserving weighted remap configured via YAML method: "weighted")
    pRaster.run_remap()
    pRaster.report_outputs()

    # 6. Visualize remapped target mesh
    pRaster.visualize_target_mesh()

    # 7. Create optional 3D animation if configured
    if config.get("visualization", {}).get("animation"):
        pRaster.visualize_target_mesh(iFlag_create_animation=True)

    # 8. Clean up temporary files
    pRaster.cleanup()
    print("done")


if __name__ == "__main__":
    main()
