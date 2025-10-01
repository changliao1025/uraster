"""
uraster - A Python package for unstructured raster processing and remapping
"""

__version__ = "0.1.0"
__author__ = "Chang Liao"
__email__ = "changliao.climate@gmail.com"

# Import main classes for easier access
from .classes.uraster import uraster
from .classes.sraster import sraster

__all__ = ["uraster", "sraster"]