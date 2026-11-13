import os
import tempfile
from pathlib import Path

import pytest

from uraster import (
    DEFAULT_CONFIG,
    URasterWorkflow,
    URasterWorkflowSpec,
    compile_remap_workflow,
    create_workflow,
    dump_workflow_config,
    load_workflow_config,
)


def test_workflow_spec_roundtrip():
    data = {
        "workflow": {"name": "test_run", "debug": True},
        "mesh": {"source_mesh": "mesh.geojson", "cellid_field": "cellid"},
        "rasters": {"source_rasters": ["r1.tif", "r2.tif"], "is_discrete": True},
        "remap": {"method": "weighted", "global": True, "polar": False},
        "outputs": {"output_workspace": "./out", "target_mesh": "res.geojson"},
        "visualization": {"enabled": False},
    }
    spec = URasterWorkflowSpec.from_dict(data)
    assert spec.workflow["name"] == "test_run"
    assert spec.rasters["is_discrete"] is True

    exported = spec.to_dict()
    assert exported == data


def test_compile_remap_workflow():
    data = {
        "mesh": {"source_mesh": "/path/to/mesh.geojson", "cellid_field": "fid"},
        "rasters": {"source_rasters": ["/path/to/r.tif"], "is_discrete": True},
        "remap": {"method": "weighted", "global": True, "polar": False},
        "outputs": {
            "output_workspace": "/path/to/out",
            "target_mesh": "remapped.geojson",
        },
    }
    compiled = compile_remap_workflow(data)
    assert compiled["sFilename_source_mesh"] == "/path/to/mesh.geojson"
    assert compiled["sField_unique_id"] == "fid"
    assert compiled["aFilename_source_raster"] == ["/path/to/r.tif"]
    assert compiled["sFilename_target_mesh"] == os.path.join(
        "/path/to/out", "remapped.geojson"
    )
    assert compiled["iFlag_remap_method"] == 3
    assert compiled["iFlag_discrete"] == 1
    assert compiled["iFlag_global"] == 1
    assert compiled["iFlag_polar"] == 0


def test_load_and_dump_workflow_config_with_comments(tmp_path):
    yaml_content = """
    # ==============================================================================
    # Test YAML workflow with comments
    # ==============================================================================
    workflow:
      name: my_test_case   # Case identifier
      debug: true          # Verbose output

    mesh:
      source_mesh: mesh.geojson # Relative path
      cellid_field: cellid

    rasters:
      source_rasters:
        - elevation.tif     # Source DEM tile
      is_discrete: false

    remap:
      method: nearest       # 1=nearest, 2=bilinear, 3=weighted
      global: true

    outputs:
      output_workspace: ./results
      target_mesh: output.geojson
    """
    config_file = tmp_path / "workflow_config.yaml"
    config_file.write_text(yaml_content, encoding="utf-8")

    # Test load with relative path resolution
    cfg = load_workflow_config(config_file)
    assert cfg["workflow"]["name"] == "my_test_case"
    assert cfg["workflow"]["debug"] is True
    # Verify relative paths were resolved against config_file.parent
    assert Path(cfg["mesh"]["source_mesh"]) == (tmp_path / "mesh.geojson").resolve()
    assert (
        Path(cfg["rasters"]["source_rasters"][0])
        == (tmp_path / "elevation.tif").resolve()
    )
    assert Path(cfg["outputs"]["output_workspace"]) == (tmp_path / "results").resolve()

    # Test directory auto-discovery
    cfg_dir = load_workflow_config(tmp_path)
    assert cfg_dir["workflow"]["name"] == "my_test_case"

    # Test dumping to JSON
    json_out = tmp_path / "exported.json"
    dump_workflow_config(cfg, json_out)
    assert json_out.is_file()

    loaded_json = load_workflow_config(json_out)
    assert loaded_json["workflow"]["name"] == "my_test_case"


def test_create_workflow_factory(tmp_path):
    yaml_content = """
    mesh:
      source_mesh: input_mesh.geojson
    rasters:
      source_rasters: [raster.tif]
    remap:
      method: weighted
    outputs:
      output_workspace: ./run_out
      target_mesh: out.geojson
    """
    config_file = tmp_path / "workflow_config.yaml"
    config_file.write_text(yaml_content, encoding="utf-8")

    wf = create_workflow(config_file)
    assert isinstance(wf, URasterWorkflow)
    assert wf.compiled_config["iFlag_remap_method"] == 3
    assert wf.output_workspace == str((tmp_path / "run_out").resolve())
