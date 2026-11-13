# Two-Tier Compiler Pattern for URaster Workflows

> **Status:** ACTIVE SPECIFICATION & IMPLEMENTATION
> - **Tier 1 (Front-End Specification):** YAML / JSON workflow specification via [`URasterWorkflowSpec`](uraster/config/workflow_spec.py) and [`load_workflow_config()`](uraster/config/config_manager.py).
> - **Compiler Boundary:** [`compile_remap_workflow()`](uraster/config/config_manager.py) compiles high-level human domain intent into backend execution parameters (`aConfig`).
> - **Tier 2 (Execution Orchestration):** [`URasterWorkflow`](uraster/config/workflow.py) coordinates execution lifecycle (`setup()`, `run()`, `visualize()`) over the core [`URaster`](uraster/classes/uraster.py) engine.

---

## 1. Motivation

Scientific raster-to-mesh remapping requires balancing two needs:
1. **User Friendliness & Maintainability:** Users need human-readable configuration files with inline comments (`# ...`) to document data origins, coordinate choices, and method tradeoffs.
2. **Execution Precision:** The underlying geospatial processing engines require precise numerical flags, array bounds, and formatted driver strings.

Following the Two-Tier Compiler Pattern designed in `pyearthmesh`, `uraster` separates high-level workflow intent from low-level runtime configuration.

---

## 2. Architecture Overview

```text
User / Configuration File (workflow_config.yaml with comments)
        │
        ▼
Tier 1: load_workflow_config() / URasterWorkflowSpec
  ├── Discovers workflow_config.yaml / .yml / .json
  ├── Resolves relative paths against the YAML directory
  ├── Resolves OS-specific paths (Linux / Darwin / Windows)
  └── Merges defaults cleanly
        │
        ▼
Compiler Boundary: compile_remap_workflow()
  ├── Maps declarative methods ('weighted' -> 3, 'nearest' -> 1)
  ├── Maps discrete data flags (is_discrete: true -> 1)
  ├── Joins output workspaces and mesh targets
  └── Emits backend-native aConfig dictionary
        │
        ▼
Tier 2: URasterWorkflow
  ├── workflow.setup()       -> validates inputs, ensures output workspace
  ├── workflow.run()         -> executes URaster.run_remap()
  └── workflow.visualize()   -> renders static images and animations
```

---

## 3. Tier 1 YAML Specification Reference

```yaml
# ==============================================================================
# URaster Workflow Specification
# ==============================================================================
workflow:
  name: "global_dem_remap"           # Case identifier
  date: null                         # Defaults to today (YYYYMMDD) if null
  debug: false                       # Diagnostic logging

mesh:
  source_mesh: "mesh.geojson"        # Source mesh file (GeoJSON, FlatGeobuf, etc.)
  cellid_field: "cellid"             # Unique ID field name in mesh
  mesh_type: "generic"               # Mesh family hint (mpas, dggs, rhealpix, tin)

rasters:
  source_rasters:
    - "elevation_tile1.tif"          # One or more source raster tiles
  is_discrete: false                 # false for continuous, true for categorical

remap:
  method: "weighted"                 # 'nearest', 'bilinear', or 'weighted'
  global: true                       # Global domain coverage
  polar: false                       # Polar mesh handling
  resolution_ratio_threshold: 3.0    # Ratio to trigger weighted averaging

outputs:
  output_workspace: "results"        # Output directory
  target_mesh: "uraster.geojson"     # Remapped vector file name

visualization:
  enabled: true                      # Whether to render output visualizations
  source_mesh_image: "mesh.jpg"
  target_mesh_image: "uraster.png"
  animation: "uraster.gif"
  focus_coordinates: null            # Optional [lon, lat] center
```

---

## 4. Python API Usage

### Quick Execution
```python
import uraster

# Automatically discovers workflow_config.yaml in the directory and runs
workflow = uraster.create_workflow("path/to/workflow_dir")
workflow.run()
```

### Inspecting Compiled Config
```python
import uraster

spec = uraster.load_workflow_config("workflow_config.yaml")
backend_config = uraster.compile_remap_workflow(spec)
print(backend_config["iFlag_remap_method"])  # 3 (weighted)
```
