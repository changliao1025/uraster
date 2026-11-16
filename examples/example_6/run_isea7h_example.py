import os

from uraster import URaster, compile_remap_workflow, load_workflow_config


def main():
    # 1. Load Tier 1 configuration from workflow_config.yaml
    sFilename_yaml = os.path.join(os.path.dirname(__file__), "workflow_config.yaml")
    config = load_workflow_config(sFilename_yaml)

    # 2. Compile Tier 1 specification into backend execution parameters
    aConfig = compile_remap_workflow(config)
    pRaster = URaster(aConfig)

    # 3. Setup and report input metadata
    pRaster.setup()
    pRaster.report_inputs()

    # 4. Visualize source mesh
    pRaster.visualize_source_mesh()

    # 5. Run remap (polar handling configured via YAML polar: true)
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
