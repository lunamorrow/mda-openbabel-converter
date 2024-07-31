# Testing OpenBabel and Pybel

import MDAnalysis as mda
import openbabel
from openbabel import openbabel as ob
from openbabel.openbabel import OBMol, OBConversion, GetSymbol, OBMolAtomIter
from mda_openbabel_converter.data.files import CRN

import mda_openbabel_converter
import pytest  # version 8.2.2
import sys
import numpy as np
from numpy.testing import assert_equal, assert_allclose
from MDAnalysisTests.topology.base import ParserBase
from MDAnalysis.core.topology import Topology

import mda_openbabel_converter.OpenBabel

class TestOpenBabelReader(object):
    reader = mda_openbabel_converter.OpenBabel.OpenBabelReader

    # @pytest.mark.parametrize("obmol, n_frames", [
    #     (CRN, 1),
    #    # ("smiles_mol", 3),
    # ], indirect=["obmol"])

    # @pytest.fixture()
    # def obmol(self):
    #     OBConversion().SetInFormat("pdb")
    #     mol = OBMol()
    #     OBConversion.ReadFile(mol, CRN)
    #     return mol

    # @pytest.fixture()
    # def obmol(self):
    #     obConversion = ob.OBConversion()
    #     obConversion.SetInFormat("pdb")
    #     mol = OBMol()
    #     obConversion.ReadFile(mol, CRN)
    #     return mol
    
    @pytest.fixture()
    def n_frames(self):
        return 1

    def test_coordinates(self, n_frames):

        obconversion = ob.OBConversion()
        obconversion.SetInFormat("pdb")
        obmol = ob.OBMol()
        obconversion.ReadFile(obmol, "mda_openbabel_converter/data/1crn.pdb")
        universe = mda.Universe(obmol)
        assert universe.trajectory.n_frames == n_frames
        expected = np.array([
            [atom.GetVector().GetX(),
            atom.GetVector().GetY(),
            atom.GetVector().GetZ()] for atom in OBMolAtomIter(obmol)],
            dtype=np.float32)
        assert_equal(expected, universe.trajectory.coordinate_array[0])
