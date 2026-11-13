"""Workflow specification module for uraster.

Captures high-level, human-editable domain intent (Tier 1) before compilation
into backend execution parameters (Tier 2).
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


@dataclass
class URasterWorkflowSpec:
    """High-level domain specification for a raster-to-mesh remapping workflow.

    Represents Tier 1 in the two-tier compiler architecture. It captures
    user intent cleanly and agnostically without hardcoding low-level
    internal runtime keys.
    """

    workflow: Dict[str, Any] = field(default_factory=dict)
    mesh: Dict[str, Any] = field(default_factory=dict)
    rasters: Dict[str, Any] = field(default_factory=dict)
    remap: Dict[str, Any] = field(default_factory=dict)
    outputs: Dict[str, Any] = field(default_factory=dict)
    visualization: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "URasterWorkflowSpec":
        """Instantiate specification from a parsed dictionary."""
        if not isinstance(data, dict):
            raise ValueError("Workflow specification must be a dictionary")

        return cls(
            workflow=dict(data.get("workflow") or {}),
            mesh=dict(data.get("mesh") or {}),
            rasters=dict(data.get("rasters") or {}),
            remap=dict(data.get("remap") or {}),
            outputs=dict(data.get("outputs") or {}),
            visualization=dict(data.get("visualization") or {}),
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert specification back to a clean sectioned dictionary."""
        return {
            "workflow": self.workflow,
            "mesh": self.mesh,
            "rasters": self.rasters,
            "remap": self.remap,
            "outputs": self.outputs,
            "visualization": self.visualization,
        }
