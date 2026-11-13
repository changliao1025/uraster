"""
uraster - A Python package for unstructured raster processing and remapping
"""

__version__ = "0.1.9"
__author__ = "Chang Liao"
__email__ = "changliao.climate@gmail.com"

from .classes.sraster import SRaster

# Import main classes for easier access
from .classes.uraster import URaster

# Two-tier workflow configuration and orchestration API
from .config import (
    DEFAULT_CONFIG,
    URasterWorkflow,
    URasterWorkflowSpec,
    compile_remap_workflow,
    create_workflow,
    dump_workflow_config,
    load_workflow_config,
)

__all__ = [
    "URaster",
    "SRaster",
    "URasterWorkflow",
    "URasterWorkflowSpec",
    "create_workflow",
    "load_workflow_config",
    "dump_workflow_config",
    "compile_remap_workflow",
    "DEFAULT_CONFIG",
]
