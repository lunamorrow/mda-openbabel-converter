"""
mda-openbabel-converter
A package to convert between MDAnalysis and OpenBabel Objects
"""

# Add imports here
from importlib.metadata import version
from .OpenBabel import OpenBabelReader
from .OpenBabel import OpenBabelConverter
# from .OpenBabelParser import OpenBabelTopologyParser

__version__ = version("mda_openbabel_converter")
