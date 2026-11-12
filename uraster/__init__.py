"""
uraster - A Python package for unstructured raster processing and remapping
"""

__version__ = "0.1.9"
__author__ = "Chang Liao"
__email__ = "changliao.climate@gmail.com"

from .classes.sraster import SRaster

# Import main classes for easier access
from .classes.uraster import URaster

__all__ = ["URaster", "SRaster"]
