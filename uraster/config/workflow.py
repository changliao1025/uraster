"""Workflow orchestration adapter for uraster.

Tier 2 execution wrapper that compiles specifications and coordinates
URaster lifecycle methods: setup(), run_remap(), and visualization.
"""

import os
from pathlib import Path
from typing import Any, Dict, Optional, Union

from uraster.classes.uraster import URaster
from uraster.config.config_manager import (
    compile_remap_workflow,
    load_workflow_config,
)
from uraster.config.workflow_spec import URasterWorkflowSpec
from uraster.utility import setup_logger

logger = setup_logger(__name__.split(".")[-1])


class URasterWorkflow:
    """High-level workflow coordinator for an unstructured raster remap job."""

    def __init__(
        self,
        config_or_spec: Union[str, Path, Dict[str, Any], URasterWorkflowSpec],
    ):
        if isinstance(config_or_spec, (str, Path)):
            self.raw_config = load_workflow_config(config_or_spec)
        elif isinstance(config_or_spec, URasterWorkflowSpec):
            self.raw_config = config_or_spec.to_dict()
        elif isinstance(config_or_spec, dict):
            self.raw_config = dict(config_or_spec)
        else:
            raise TypeError(
                "config_or_spec must be a path (str/Path), URasterWorkflowSpec, or dict"
            )

        self.compiled_config = compile_remap_workflow(self.raw_config)
        self.output_workspace = (
            self.raw_config.get("outputs", {}).get("output_workspace")
            or os.path.dirname(self.compiled_config.get("sFilename_target_mesh", ""))
            or "./output"
        )
        self.engine = URaster(self.compiled_config)

    def setup(self, iFlag_verbose_in: bool = True) -> bool:
        """Create output workspace and validate engine inputs."""
        os.makedirs(self.output_workspace, exist_ok=True)
        res = self.engine.setup(iFlag_verbose_in=iFlag_verbose_in)
        if res and iFlag_verbose_in:
            self.engine.report_inputs()
        return res

    def run(self, iFlag_verbose_in: bool = True) -> bool:
        """Execute the remapping operation."""
        self.setup(iFlag_verbose_in=iFlag_verbose_in)
        logger.info(
            f"Starting remap workflow: {self.compiled_config.get('sFilename_source_mesh')} "
            f"-> {self.compiled_config.get('sFilename_target_mesh')}"
        )
        # Check if discrete or weighted flags should be passed explicitly
        is_discrete = bool(self.compiled_config.get("iFlag_discrete", 0))
        is_weighted = self.compiled_config.get("iFlag_remap_method", 1) == 3

        self.engine.run_remap(
            iFlag_discrete_in=is_discrete,
            iFlag_weighted_average_in=is_weighted,
            iFlag_verbose_in=iFlag_verbose_in,
        )

        # Optional visualization
        vis_config = self.raw_config.get("visualization", {})
        if vis_config.get("enabled", False):
            self.visualize(vis_config)

        self.engine.cleanup()
        return True

    def visualize(self, vis_config: Optional[Dict[str, Any]] = None) -> None:
        """Render optional static images and 3D animations."""
        cfg = vis_config or self.raw_config.get("visualization", {})
        focus_coords = cfg.get("focus_coordinates")
        lon_focus = focus_coords[0] if focus_coords else None
        lat_focus = focus_coords[1] if focus_coords else None

        # Source mesh
        mesh_img = cfg.get("source_mesh_image")
        if mesh_img:
            out_mesh = (
                mesh_img
                if os.path.isabs(mesh_img)
                else os.path.join(self.output_workspace, mesh_img)
            )
            try:
                self.engine.visualize_source_mesh(
                    sFilename_out=out_mesh,
                    dLongitude_focus_in=lon_focus,
                    dLatitude_focus_in=lat_focus,
                )
            except Exception as e:
                logger.warning(f"Could not render source mesh visualization: {e}")

        # Target mesh
        target_img = cfg.get("target_mesh_image")
        colormap = cfg.get("colormap", "terrain")
        if target_img:
            out_target = (
                target_img
                if os.path.isabs(target_img)
                else os.path.join(self.output_workspace, target_img)
            )
            try:
                self.engine.visualize_target_mesh(
                    sFilename_out=out_target,
                    sColormap=colormap,
                    dLongitude_focus_in=lon_focus,
                    dLatitude_focus_in=lat_focus,
                )
            except Exception as e:
                logger.warning(f"Could not render target mesh visualization: {e}")

        # Animation
        anim_file = cfg.get("animation")
        if anim_file:
            out_anim = (
                anim_file
                if os.path.isabs(anim_file)
                else os.path.join(self.output_workspace, anim_file)
            )
            try:
                anim_format = "mp4" if out_anim.endswith(".mp4") else "gif"
                self.engine.visualize_target_mesh(
                    sFilename_out=out_anim,
                    sColormap=colormap,
                    dLongitude_focus_in=lon_focus,
                    dLatitude_focus_in=lat_focus,
                    iFlag_create_animation=True,
                    iAnimation_frames=360,
                    sAnimation_format=anim_format,
                )
            except Exception as e:
                logger.warning(f"Could not render animation visualization: {e}")


def create_workflow(
    config_or_path: Union[str, Path, Dict[str, Any], URasterWorkflowSpec],
) -> URasterWorkflow:
    """Factory function to build a URasterWorkflow from a file, dict, or spec."""
    return URasterWorkflow(config_or_path)
