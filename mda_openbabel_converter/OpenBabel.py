"""
Documentation...
"""

import MDAnalysis as mda
from MDAnalysis.converters.base import ConverterBase
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.core.groups import AtomGroup

try:
    from openbabel import openbabel as OB
    from openbabel import OBMol
except ImportError:
    print("Cannot find openbabel, install with 'pip install openbabel==2.4.0'")

class OpenBabelReader(MemoryReader):
    """
    Convert an OpenBabel OBMol (from the file) to a MDAnalysis AtomGroup
    """
    @staticmethod
    def _format_hint(thing):
        """
        Base function to check if the reader can actually read this “thing” 
        (i.e., is it a file that can be converted to an OpenBabel OBMol?)
        """
        try:
            import openbabel as OB
        except ImportError:
            return False
        else:
            return isinstance(thing, OB.OBMol)

    def __init__(self, filename: OBMol, **kwargs):
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