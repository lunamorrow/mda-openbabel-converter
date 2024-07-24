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
from enum import StrEnum

HAS_OBABEL = False
NEUTRON_MASS = 1.008

try:
    from openbabel import openbabel as ob
    from openbabel import OBMol
    from openbabel import OBElementTable
    HAS_OBABEL = True
except ImportError:
    # import breaks with version 3.x
    warnings.warn("Cannot find openbabel, install with `mamba install -c "
                  "conda-forge openbabel==2.4.0`")


class StereoEnum(StrEnum):
    POSITIVE = "+"
    NEGATIVE = "-"
    NONE = ""


class OpenBabelParser(TopologyReaderBase):
    """
    Inherits from TopologyReaderBase and converts an OpenBabel OBMol to a
    MDAnalysis Topology or adds it to a pre-existing Topology. This parser
    does not work in the reverse direction.
    """
    format = 'OPENBABEL'

    @staticmethod
    def _format_hint(thing):
        """
        Base function to check if the parser can actually parse this “thing”
        (i.e., is it a valid OpenBabel OBMol that can be converted to a
        MDAnalysis Topology?)
        """
        if HAS_OBABEL is False:
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
        mol = self.filename

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
        chainids = []
        icodes = []

        if mol.Empty():
            return Topology(n_atoms=0,
                            n_res=0,
                            n_seg=0,
                            attrs=None,
                            atom_resindex=None,
                            residue_segindex=None)

        for atom in ob.OBMolAtomIter(mol):
            # Atom name set with element and id, as name not supported by OpenBabel
            name = "%s%d" % (OBElementTable().GetSymbol(atom.GetAtomicNum()),
                             atom.GetIdx())
            names.append(name)
            atomtypes.append(atom.GetType())
            ids.append(atom.GetIdx())
            masses.append(atom.GetExactMass())
            if abs(atom.GetExactMass()-atom.GetAtomicMass()) >= NEUTRON_MASS:
                warnings.warn(
                    f"Exact mass and atomic mass of atom ID: {atom.GetIdx()}"
                    " are more than 1.008 AMU different. Be aware of isotopes,"
                    " which are NOT flagged by MDAnalysis.")
            charges.append(atom.GetPartialCharge())

            # convert atomic number to element
            elements.append(OBElementTable().GetSymbol(atom.GetAtomicNum()))

            # only for PBD and MOL2
            if atom.HasResidue():
                resid = atom.GetResidue()
                resnums.append(resid.GetNum())
                resnames.append(resid.GetName())
                chainids.append(resid.GetChain())
                icodes.append(resid.GetInsertionCode())

            chirality = StereoEnum.NONE
            if atom.IsPositiveStereo():
                chirality = StereoEnum.POSITIVE
            elif atom.IsNegativeStereo():
                chirality = StereoEnum.NEGATIVE
            chiralities.append(chirality)

            aromatics.append(atom.IsAromatic())

        # make Topology attributes
        attrs = []
        n_atoms = len(ids)

        if resnums and (len(resnums) != len(ids)):
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
        bond_orders = []
        for bond_idx in range(0, mol.NumBonds()):
            bond = mol.GetBond(bond_idx)
            bonds.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
            bond_orders.append(float(bond.GetBondOrder()))
        attrs.append(Bonds(bonds, order=bond_orders))

        # * Optional attributes *
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
            pass  # no guesser yet

        # Residue
        if resnums:
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
