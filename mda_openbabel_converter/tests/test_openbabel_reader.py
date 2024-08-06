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

    @pytest.fixture()
    def n_frames(self):
        return 1

    def test_coordinates(self, n_frames):
        obconversion = ob.OBConversion()
        obconversion.SetInFormat("pdb")
        obmol = ob.OBMol()
        obconversion.ReadFile(obmol, CRN.as_posix())
        universe = mda.Universe(obmol)
        assert universe.trajectory.n_frames == n_frames
        expected = np.array([
            [(coords := atom.GetVector()).GetX(),
            coords.GetY(),
            coords.GetZ()] for atom in OBMolAtomIter(obmol)],
            dtype=np.float32)
        assert_equal(expected, universe.trajectory.coordinate_array[0])

    def test_no_coordinates(self):
        obConversion = ob.OBConversion()
        obConversion.SetInFormat("smi")
        mol = OBMol()
        obConversion.ReadString(mol, "C1=CC=CS1")
        mol.AddHydrogens()
        universe = mda.Universe(mol)
        n_atoms = mol.NumAtoms()
        expected = coordinates = np.empty((1, n_atoms, 3), dtype=np.float32)
        coordinates[:] = np.nan
        assert_equal(universe.trajectory.coordinate_array, expected)
