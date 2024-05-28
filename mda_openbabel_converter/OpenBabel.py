"""
Documentation...
"""

import MDAnalysis as mda
from MDAnalysis.converters.base import ConverterBase
from MDAnalysis.coordinates.base import SingleFrameReaderBase
from MDAnalysis.core.groups import AtomGroup

try:
    import openbabel as OB
    from openbabel import OBMol
except ImportError:
    print("Cannot find openbabel, install with 'pip install openbabel==2.4.0'")

class OpenBabelReader(SingleFrameReaderBase):
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
        self.atoms = []
        self.n_atoms = 0
        self.residues = []
        self.n_residues = 0
        self.segments = []
        self.n_segments = 0

        obmol = filename

        # Atoms
        names = [] 
        chiralities = []
        resnums = []
        resnames = []
        elements = []
        masses = []
        charges = []
        aromatics = []
        ids = []
        atomtypes = []
        segids = []
        altlocs = []
        chainids = []
        icodes = []
        occupancies = []
        tempfactors = []

        for i in range(1, obmol):
            atom = obmol.GetAtomById(i-1)
            # need to add handling incase attributes are invalid or null in OBMol
            names.append(atom.GetType()) #char


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