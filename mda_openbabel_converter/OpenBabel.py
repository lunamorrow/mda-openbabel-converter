"""
Documentation...
"""

import MDAnalysis as mda
from MDAnalysis.converters.base import ConverterBase
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.core.groups import AtomGroup
import numpy as np
import warnings

HAS_OBABEL = False

try:
    from openbabel import openbabel as ob
    from openbabel.openbabel import OBMol, OBAtom, OBMolAtomIter
    HAS_OBABEL = True
except ImportError:
    warnings.warn("Cannot find openbabel, install with `mamba install -c "
                  "conda-forge openbabel`")


class OpenBabelReader(MemoryReader):
    """
    Inherits from MemoryReader and converts OpenBabel OBMol Coordinates to a
    MDAnalysis Trajectory which is used to build a Universe. This reader
    does not work in the reverse direction.
    """
    format = 'OPENBABEL'

    # Structure.coordinates always in Angstrom
    units = {'time': None, 'length': 'Angstrom'}

    @staticmethod
    def _format_hint(thing):
        """
        Base function to check if the reader can actually read this “thing”
        (i.e., is it a file that can be converted to an OpenBabel OBMol?)
        """
        if HAS_OBABEL is False:
            return False
        else:
            return isinstance(thing, OBMol)

    def __init__(self, filename: OBMol, **kwargs):
        """
        Converts file to OBMol to AtomGroup
        """
        n_atoms = filename.NumAtoms()
        coordinates = np.array([
            [(coords := atom.GetVector()).GetX(),
            coords.GetY(),
            coords.GetZ()] for atom in OBMolAtomIter(filename)],
            dtype=np.float32)
        if not np.any(coordinates):
            warnings.warn("No coordinates found in the OBMol")
            coordinates = np.empty((1, n_atoms, 3), dtype=np.float32)
            coordinates[:] = np.nan
        super(OpenBabelReader, self).__init__(coordinates, order='fac', **kwargs)

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