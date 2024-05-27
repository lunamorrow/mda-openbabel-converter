"""
Documentation...
"""

import MDAnalysis as mda
from MDAnalysis.converters.base import ConverterBase
from MDAnalysis.coordinates.base import SingleFrameReaderBase

class OpenBabelReader(SingleFrameReaderBase):
    """
    Convert an OpenBabel OBMol (from the file) to a MDAnalysis AtomGroup
    """
    def _format_hint(thing):
        """
        Base function to check if the reader can actually read this “thing” 
        (i.e., is it a file that can be converted to an OpenBabel OBMol?)
        """
        pass

    def __init__(self, filename, **kwargs):
        """
        Converts file to OBMol to AtomGroup
        """
        pass

class OpenBabelConverter(ConverterBase):
    """
    Convert a MDAnalysis AtomGroup to an OpenBabel OBMol
    """
    def __repr__(self, **kwargs):
        """
        String representation of the object (defined in Base Class)
        """
        pass

    def convert(self, obj, **kwargs):
        """
        Converts AtomGroup to OBMol
        """
        pass

    # add getter and setter methods as required

    # add helper methods