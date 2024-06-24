"""
Documentation...
"""

import MDAnalysis as mda
from MDAnalysis.topology.base import TopologyReaderBase
from MDAnalysis.core.topology import Topology
from MDAnalysis.converters.base import ConverterBase
import warnings

try:
    import openbabel as OB
    from openbabel import OBMol
    from openbabel import OBElementTable
except ImportError:
    print("Cannot find openbabel, install with 'pip install openbabel==2.4.0'")

class OpenBabelParser(TopologyReaderBase):
    """
    Inherits from TopologyReaderBase and converts an OpenBabel OBMol to a 
    MDAnalysis Topology or adds it to a pre-existing Topology. This parser will 
    does not work in the reverse direction.
    """

    @staticmethod
    def _format_hint(thing):
        """
        Base function to check if the parser can actually parse this “thing”
        (i.e., is it a valid OpenBabel OBMol that can be converted to a
        MDAnalysis Topology?)
        """
        try:
            import openbabel as OB
        except ImportError:
            return False
        else:
            return isinstance(thing, OB.OBMol)

    def parse(self, **kwargs):
        """
        Accepts an OpenBabel OBMol and returns a MDAnalysis Topology. Will need
        to extract the number of atoms, number of residues, number of segments,
        atom_residue index, residue_segment index and all of the atom's
        relevant attributes from the OBMol to initialise a new Topology.
        """
        format = 'OPENBABEL'

        self.atoms = []
        self.n_atoms = 0
        self.residues = []
        self.n_residues = 0
        self.segments = []
        self.n_segments = 0

        obmol = self.filename

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
        tempfactors = []  # B factor; not supported by OB

        if obmol.GetFirstAtom().equals(None):
            return Topology(n_atoms=0,
                            n_res=0,
                            n_seg=0,
                            attrs=None,
                            atom_resindex=None,
                            residue_segindex=None)

        for atom in OB.OBMolAtomIter(obmol):
            # need to add handling incase attributes are invalid or null in OBMol
            # names.append(atom.GetType()) #char -> nothing for name in OBMol? Is name required to make MDA Atom?
            atomtypes.append(atom.GetType())  # char
            ids.append(atom.GetIdx()) #int
            masses.append(atom.GetExactMass())  # double -> what about atom.GetAtomicMass()??; which is better?
            if not atom.GetExactMass().equals(atom.GetAtomicMass()):
                warnings.warn(
                    f"Exact mass and atomic mass of atom (ID: {atom.GetIdx})
                    not equal. Be aware of isotopes, which are NOT supported
                    by MDAnalysis.")
            charges.append(atom.GetPartialCharge()) #int (or use atom.GetFormalCharge()?)

            # convert atomic number to element
            elements.append(OBElementTable.GetSymbol(atom.GetAtomicNumber())) #char

            resid = atom.GetResidue() # null if no residue
            if resid.equals(None):
                warnings.warn(
                    f"No residue is defined for atom (ID: {atom.GetIdx}). 
                    Please set with 'SETTTING METHOD FOR OBMOL'" # TO DO
                )
                # if residue null, will need to assign w MDAnalysis...
            else:
                resnums.append(resid.GetNum()) # TO DO: check if start at 0 or 1
                resnames.append(resid.GetName())

            # don't need to check null case, as know assigned to OBMol we're currently parsing
            # but, NEED TO HANDLE ADDING MULTIPLE SEGIDS/OBMOLS WHEN CONVERTING TO UNIVERSE AND ADDING TOGETHER... (should be ok, check w tests)
            segids.append(atom.GetParent())

            chiralities.append(atom.IsChiral()) #boolean
            aromatics.append(atom.IsAromatic()) #boolean

        pass