"""Workflow configuration manager and two-tier compiler for uraster.

Provides functions to load, validate, dump, and compile YAML/JSON workflow
configurations into backend execution parameters for the URaster engine.
"""

import datetime
import json
import os
import platform
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import yaml

from uraster.config.workflow_spec import URasterWorkflowSpec

DEFAULT_CONFIG: Dict[str, Dict[str, Any]] = {
    "workflow": {
        "name": "uraster_workflow",
        "date": None,
        "debug": False,
        "example_number": None,  # If set, uses Pooch to download example dataset
    },
    "mesh": {
        "source_mesh": "",
        "cellid_field": "cellid",
        "mesh_type": "generic",
    },
    "rasters": {
        "source_rasters": [],
        "is_discrete": False,
    },
    "remap": {
        "method": "nearest",  # 'nearest', 'bilinear', or 'weighted'
        "global": True,
        "polar": False,
        "resolution_ratio_threshold": 3.0,
    },
    "outputs": {
        "output_workspace": "./output",
        "target_mesh": "uraster.geojson",
    },
    "visualization": {
        "enabled": False,
        "source_mesh_image": "mesh.jpg",
        "target_mesh_image": "uraster.png",
        "animation": "uraster.gif",
        "focus_coordinates": None,
    },
}

METHOD_MAP: Dict[str, int] = {
    "nearest": 1,
    "nearest_neighbor": 1,
    "bilinear": 2,
    "weighted": 3,
    "weighted_average": 3,
    "weighted_area": 3,
}


def _deep_merge(target: Dict[str, Any], source: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively merge source dictionary into a deepcopy of target."""
    result = deepcopy(target)
    for key, val in source.items():
        if key in result and isinstance(result[key], dict) and isinstance(val, dict):
            result[key] = _deep_merge(result[key], val)
        else:
            result[key] = deepcopy(val)
    return result


def _resolve_platform_value(val: Any) -> Any:
    """Resolve OS-specific dictionary (e.g. {'Linux': '...', 'Darwin': '...', 'default': '...'})."""
    if isinstance(val, dict):
        current_os = platform.system()
        if current_os in val:
            return val[current_os]
        if "default" in val:
            return val["default"]
    return val


def _resolve_paths(data: Any, base_dir: Path) -> Any:
    """Recursively resolve platform values and relative paths against base_dir."""
    if isinstance(data, dict):
        resolved: Dict[str, Any] = {}
        for k, v in data.items():
            v_plat = _resolve_platform_value(v)
            if (
                isinstance(v_plat, str)
                and v_plat.strip()
                and (
                    "mesh" in k
                    or "raster" in k
                    or "workspace" in k
                    or "file" in k
                    or "path" in k
                )
            ):
                # Path resolution
                path_obj = Path(v_plat)
                if not path_obj.is_absolute():
                    resolved[k] = str((base_dir / path_obj).resolve())
                else:
                    resolved[k] = str(path_obj)
            elif isinstance(v_plat, list):
                resolved[k] = [
                    (
                        str((base_dir / Path(item)).resolve())
                        if isinstance(item, str)
                        and item.strip()
                        and not Path(item).is_absolute()
                        else item
                    )
                    for item in v_plat
                ]
            else:
                resolved[k] = _resolve_paths(v_plat, base_dir)
        return resolved
    return data


def find_workflow_config_file(path_or_dir: Union[str, Path]) -> Path:
    """Locate workflow configuration file from a directory, exact file, or stem."""
    p = Path(path_or_dir).expanduser().resolve()
    if p.is_dir():
        for candidate in [
            "workflow_config.yaml",
            "workflow_config.yml",
            "workflow_config.json",
        ]:
            c = p / candidate
            if c.is_file():
                return c
        raise FileNotFoundError(f"No workflow_config (.yaml, .yml, .json) found in {p}")

    if p.is_file():
        return p

    # Check potential extensions if passed without extension
    for ext in [".yaml", ".yml", ".json"]:
        cand = p.with_suffix(ext)
        if cand.is_file():
            return cand

    raise FileNotFoundError(f"Configuration file not found: {path_or_dir}")


def load_workflow_config(path_or_dir: Union[str, Path]) -> Dict[str, Any]:
    """Load and validate a YAML or JSON workflow configuration with comments."""
    config_file = find_workflow_config_file(path_or_dir)
    base_dir = config_file.parent

    with open(config_file, encoding="utf-8") as f:
        text = f.read()

    parsed = yaml.safe_load(text) or {}
    if not isinstance(parsed, dict):
        raise ValueError(f"Invalid workflow configuration in {config_file}")

    # Check if this is already an old-style flat aConfig dictionary
    is_flat_config = any(
        k.startswith(("sFilename_", "iFlag_", "aFilename_")) for k in parsed
    )
    if is_flat_config:
        merged = parsed
    else:
        merged = _deep_merge(DEFAULT_CONFIG, parsed)

        # If workflow has an example_number specified and Pooch is available,
        # resolve input files relative to the downloaded Pooch cache directory
        example_num = merged.get("workflow", {}).get("example_number")
        if example_num is not None:
            try:
                from uraster.utility import get_example_paths

                pooch_paths = get_example_paths(int(example_num))
                pooch_input_dir = Path(pooch_paths["input"])
                # Resolve mesh path if relative
                mesh_p = merged.get("mesh", {}).get("source_mesh", "")
                if mesh_p and not Path(mesh_p).is_absolute():
                    merged["mesh"]["source_mesh"] = str(pooch_input_dir / mesh_p)
                # Resolve raster paths if relative
                rasters_list = merged.get("rasters", {}).get("source_rasters", [])
                if isinstance(rasters_list, list):
                    merged["rasters"]["source_rasters"] = [
                        (
                            str(pooch_input_dir / r)
                            if isinstance(r, str) and not Path(r).is_absolute()
                            else r
                        )
                        for r in rasters_list
                    ]
            except Exception:
                pass

        merged = _resolve_paths(merged, base_dir)

    # Normalise date if present
    if isinstance(merged.get("workflow"), dict):
        date_val = merged["workflow"].get("date")
        if not date_val:
            merged["workflow"]["date"] = datetime.datetime.now().strftime("%Y%m%d")

    merged["_config_file"] = str(config_file)
    return merged


def dump_workflow_config(
    config: Dict[str, Any],
    file_path: Union[str, Path],
    fmt: Optional[str] = None,
) -> None:
    """Dump workflow configuration to YAML or JSON."""
    clean_config = {k: v for k, v in config.items() if not k.startswith("_")}
    path = Path(file_path).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)

    target_fmt = fmt.lower() if fmt else path.suffix.lstrip(".").lower()
    if target_fmt in ["yaml", "yml", ""]:
        with open(path, "w", encoding="utf-8") as f:
            yaml.dump(clean_config, f, default_flow_style=False, sort_keys=False)
    elif target_fmt == "json":
        with open(path, "w", encoding="utf-8") as f:
            json.dump(clean_config, f, indent=2)
    else:
        raise ValueError(f"Unsupported configuration format: {target_fmt}")


def compile_remap_workflow(
    spec_or_config: Union[URasterWorkflowSpec, Dict[str, Any]],
) -> Dict[str, Any]:
    """Compile Tier 1 workflow specification into backend-native aConfig dictionary.

    Translates human-readable domain keys (methods, discrete flag, outputs)
    into the exact dictionary format expected by the core URaster class.
    """
    if isinstance(spec_or_config, URasterWorkflowSpec):
        data = spec_or_config.to_dict()
    elif isinstance(spec_or_config, dict):
        data = spec_or_config
    else:
        raise TypeError("spec_or_config must be URasterWorkflowSpec or dict")

    # If already a legacy flat aConfig dict, pass through cleanly
    if any(k.startswith(("sFilename_", "iFlag_", "aFilename_")) for k in data):
        return dict(data)

    mesh_sec = data.get("mesh", {})
    rasters_sec = data.get("rasters", {})
    remap_sec = data.get("remap", {})
    outputs_sec = data.get("outputs", {})
    vis_sec = data.get("visualization", {})
    wf_sec = data.get("workflow", {})

    # Method flag translation
    method_val = remap_sec.get("method", "nearest")
    if isinstance(method_val, str):
        iFlag_remap_method = METHOD_MAP.get(method_val.strip().lower(), 1)
    else:
        iFlag_remap_method = int(method_val)

    # Boolean to int flags
    iFlag_discrete = 1 if rasters_sec.get("is_discrete", False) else 0
    iFlag_global = 1 if remap_sec.get("global", True) else 0
    iFlag_polar = 1 if remap_sec.get("polar", False) else 0

    # Path resolution for target mesh
    output_workspace = outputs_sec.get("output_workspace", "./output")
    target_mesh_name = outputs_sec.get("target_mesh", "uraster.geojson")
    if not os.path.isabs(target_mesh_name):
        sFilename_target_mesh = os.path.join(output_workspace, target_mesh_name)
    else:
        sFilename_target_mesh = target_mesh_name

    source_mesh = mesh_sec.get("source_mesh", "")
    source_rasters = rasters_sec.get("source_rasters", [])
    if isinstance(source_rasters, str):
        source_rasters = [source_rasters]

    aConfig: Dict[str, Any] = {
        "sFilename_source_mesh": source_mesh,
        "sField_unique_id": mesh_sec.get("cellid_field", "cellid"),
        "aFilename_source_raster": list(source_rasters),
        "sFilename_target_mesh": sFilename_target_mesh,
        "iFlag_remap_method": iFlag_remap_method,
        "iFlag_discrete": iFlag_discrete,
        "iFlag_global": iFlag_global,
        "iFlag_polar": iFlag_polar,
        "dResolution_ratio_threshold": float(
            remap_sec.get("resolution_ratio_threshold", 3.0)
        ),
        # Retain metadata and visualization sections for the workflow runner
        "_workflow": wf_sec,
        "_outputs": outputs_sec,
        "_visualization": vis_sec,
    }

    return aConfig
