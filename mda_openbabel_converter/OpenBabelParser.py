"""
Documentation...
"""

import MDAnalysis as mda
from MDAnalysis.topology.base import TopologyReaderBase, change_squash
from MDAnalysis.core.topology import Topology
from MDAnalysis.topology import guessers
from MDAnalysis.converters.base import ConverterBase
from MDAnalysis.core.topologyattrs import (
    Atomids,
    Atomnames,
    Atomtypes,
    Elements,
    Masses,
    Charges,
    Aromaticities,
    Bonds,
    Resids,
    Resnums,
    Resnames,
    RSChirality,
    Segids,
    AltLocs,
    ChainIDs,
    ICodes,
    Occupancies,
    Tempfactors,
)
import warnings
import numpy as np
import pdb #  for debugging
HAS_OBABEL=False

try:
    import openbabel as ob
    from openbabel import OBMol
    from openbabel import OBElementTable
    HAS_OBABEL=True
except ImportError:
    warnings.warn("Cannot find openbabel, install with `mamba install -c "
                  "conda-forge openbabel`")

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
        if HAS_OBABEL == False:
            return False
        else:
            return isinstance(thing, ob.OBMol)

    def parse(self, **kwargs):
        """
        Accepts an OpenBabel OBMol and returns a MDAnalysis Topology. Will need
        to extract the number of atoms, number of residues, number of segments,
        atom_residue index, residue_segment index and all of the atom's
        relevant attributes from the OBMol to initialise a new Topology.
        """
        format = 'OPENBABEL'
        mol = self.filename

        self.atoms = []
        self.n_atoms = 0
        self.residues = []
        self.n_residues = 0
        self.segments = []
        self.n_segments = 0

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

        if mol.Empty():
            return Topology(n_atoms=0,
                            n_res=0,
                            n_seg=0,
                            attrs=None,
                            atom_resindex=None,
                            residue_segindex=None)

        for atom in ob.OBMolAtomIter(mol):
            # need to add handling incase attributes are invalid or null in OBMol
            # names.append(atom.GetType()) #char -> nothing for name in OBMol? Is name required to make MDA Atom?
            atomtypes.append(atom.GetType())  # char
            ids.append(atom.GetIdx()) #int
            masses.append(atom.GetExactMass())  # double -> what about atom.GetAtomicMass()??; which is better?
            if not (atom.GetExactMass() == atom.GetAtomicMass()):
                warnings.warn(
                    f"Exact mass and atomic mass of atom (ID: {atom.GetIdx})"
                    "not equal. Be aware of isotopes, which are NOT supported"
                    "by MDAnalysis.")
            charges.append(atom.GetPartialCharge()) #int (or use atom.GetFormalCharge()?)

            # convert atomic number to element
            elements.append(OBElementTable().GetSymbol(atom.GetAtomicNum())) #char

            if atom.HasResidue():
                resid = atom.GetResidue() # null if no residue
                resnums.append(resid.GetNum()) # TO DO: check if start at 0 or 1
                resnames.append(resid.GetName())
                chainids.append(resid.GetChainNum()) # is this correct???
                icodes.append(resid.GetInsertionCode())
            else:
                warnings.warn(
                    f"No residue is defined for atom (ID: {atom.GetIdx})." 
                    "Please set with 'MDA SETTING METHOD' if required." # TO DO
                )
                resnums.append(None) # TO DO: is this best??
                resnames.append(None)
                chainids.append(None)
                icodes.append(None)
            
                 
            # don't need to check null case, as know assigned to OBMol we're currently parsing
            # but, NEED TO HANDLE ADDING MULTIPLE SEGIDS/OBMOLS WHEN CONVERTING TO UNIVERSE AND ADDING TOGETHER... (should be ok, check w tests)
            # segids.append(atom.GetParent())
            segids.append(0) # need better system!

            # TO DO: may need to create seperate for if SMILES input
            chirality = None
            if atom.IsPositiveStereo():
                chirality = "+"
            if atom.IsNegativeStereo():
                chirality = "-"
            chiralities.append(chirality)

            aromatics.append(atom.IsAromatic()) #boolean

            # altlocs.append()
            # occupancies.append()
            # tempfactors.append()


        # make Topology attributes
        attrs = []
        n_atoms = len(ids)

        if resnums and (len(resnums) != n_atoms): #resnums.__contains__(None):
            raise ValueError(
                "ResidueInfo is only partially available in the molecule."
            )

        # * Attributes always present *

        # Atom attributes
        for vals, Attr, dtype in (
            (ids, Atomids, np.int32),
            (elements, Elements, object),
            (masses, Masses, np.float32),
            (aromatics, Aromaticities, bool),
            (chiralities, RSChirality, 'U1'),
        ):
            attrs.append(Attr(np.array(vals, dtype=dtype)))

        # Bonds
        bonds = []
        bond_types = []
        bond_orders = []
        for bond_idx in range(0, mol.NumBonds()):
            bond = mol.GetBond(bond_idx)
            print(bond)
            bonds.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
            bond_orders.append(bond.GetBondOrder()) # is int, not double. Does this matter?

            # make these a dict instead?
            OB_BOND_TYPES = [
            bond.IsAromatic(),
            bond.IsAmide(),
            bond.IsPrimaryAmide(),
            bond.IsSecondaryAmide(),
            bond.IsEster(),
            bond.IsCarbonyl(),
            ]

            MDA_BOND_TYPES = [
                "aromatic",
                "amide",
                "primary amide",
                "secondary amide",
                "ester",
                "carbonyl",
            ]
            for index, b_type in enumerate(OB_BOND_TYPES):
                if b_type==True:
                    bond_types.append(MDA_BOND_TYPES[index])

        attrs.append(Bonds(bonds, types=bond_types, order=bond_orders))

        # * Optional attributes *

        # Atom name
        if names:
            attrs.append(Atomnames(np.array(names, dtype=object)))
        else:
            for atom in ob.OBMolAtomIter(mol):
                name = "%s%d" % (OBElementTable().GetSymbol(atom.GetAtomicNum()), atom.GetIdx())
                names.append(name)
            attrs.append(Atomnames(np.array(names, dtype=object)))

        # Atom type
        if atomtypes:
            attrs.append(Atomtypes(np.array(atomtypes, dtype=object)))
        else:
            atomtypes = guessers.guess_types(names)
            attrs.append(Atomtypes(atomtypes, guessed=True))

        # Partial charges
        if charges:
            attrs.append(Charges(np.array(charges, dtype=np.float32)))
        else:
            pass # no guesser yet

        # # PDB only
        # for vals, Attr, dtype in (
        #     (altlocs, AltLocs, object),
        #     (chainids, ChainIDs, object),
        #     (occupancies, Occupancies, np.float32),
        #     (tempfactors, Tempfactors, np.float32),
        # ):
        #     if vals:
        #         attrs.append(Attr(np.array(vals, dtype=dtype)))

        # Residue
        if any(resnums) and not any(val is None for val in resnums):
            resnums = np.array(resnums, dtype=np.int32)
            resnames = np.array(resnames, dtype=object)
            segids = np.array(segids, dtype=object)
            icodes = np.array(icodes, dtype=object)
            residx, (resnums, resnames, icodes, segids) = change_squash(
                (resnums, resnames, icodes, segids),
                (resnums, resnames, icodes, segids))
            n_residues = len(resnums)
            for vals, Attr, dtype in (
                (resnums, Resids, np.int32),
                (resnums.copy(), Resnums, np.int32),
                (resnames, Resnames, object),
                (icodes, ICodes, object),
            ):
                attrs.append(Attr(np.array(vals, dtype=dtype)))
        else:
            attrs.append(Resids(np.array([1])))
            attrs.append(Resnums(np.array([1])))
            residx = None
            n_residues = 1

        # Segment
        if any(segids) and not any(val is None for val in segids):
            segidx, (segids,) = change_squash((segids,), (segids,))
            n_segments = len(segids)
            attrs.append(Segids(segids))
        else:
            n_segments = 1
            attrs.append(Segids(np.array(['SYSTEM'], dtype=object)))
            segidx = None

        # create topology
        top = Topology(n_atoms, n_residues, n_segments,
                       attrs=attrs,
                       atom_resindex=residx,
                       residue_segindex=segidx)

        return top
