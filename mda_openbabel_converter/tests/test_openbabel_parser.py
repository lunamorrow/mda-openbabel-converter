# Testing OpenBabel and Pybel

import MDAnalysis as mda
import openbabel as ob
from openbabel import OBMol, OBConversion, OBElementTable
from pybel import readfile
# from openbabel.test.files import files

import mda_openbabel_converter
from mda_openbabel_converter import OpenBabelParser as OBParser
import pytest  # version 8.2.2
import sys
import numpy as np
from numpy.testing import assert_equal, assert_allclose
from MDAnalysisTests.topology.base import ParserBase
from MDAnalysis.core.topology import Topology

import mda_openbabel_converter.OpenBabelParser  # version 2.7.0

# *** can run with "python -m pytest" but not "pytest" (can't find 
# MDAnalysis) - need to fix this! ***


class OpenBabelParserBase(ParserBase):
    parser = mda_openbabel_converter.OpenBabelParser.OpenBabelParser

    expected_attrs = ['ids', 'names', 'elements', 'masses', 'aromaticities',
                      'resids', 'resnums', 'chiralities',
                      'segids', 'bonds',
                    ]

    expected_n_atoms = 0
    expected_n_residues = 0
    expected_n_segments = 0
    expected_n_bonds = 0

    def test_creates_universe(self, filename):
        u = mda.Universe(filename, format='OPENBABEL')
        assert isinstance(u, mda.Universe)

    def test_bonds_total_counts(self, top):
        if hasattr(top, 'bonds'):
            assert len(top.bonds.values) == self.expected_n_bonds


class TestOpenBabelParserEmpty(OpenBabelParserBase):
    @pytest.fixture()
    def filename(self):
        return OBMol()

    expected_attrs = []
    mandatory_attrs = []  # as not instantiated during empty Topology creation

    def test_mandatory_attributes(self, top):
        for attr in self.mandatory_attrs:
            assert (hasattr(top, attr), 
                    'Missing required attribute: {}'.format(attr))

    def test_attrs_total_counts(self, top):
        ag = mda.Universe(top).select_atoms("all")
        res = ag.residues
        seg = ag.segments
        assert len(ag) == self.expected_n_atoms
        assert len(res) == self.expected_n_residues
        assert len(seg) == self.expected_n_segments


class TestOpenBabelParserSMILES(OpenBabelParserBase):
    expected_attrs = OpenBabelParserBase.expected_attrs + ['charges']

    @pytest.fixture()
    def filename(self):
        obConversion = ob.OBConversion()
        obConversion.SetInFormat("smi")
        mol = OBMol()
        obConversion.ReadString(mol, "C1=CC=CS1")
        mol.AddHydrogens()
        return mol

    @pytest.fixture()
    def top(self, filename):
        yield self.parser(filename).parse()

    expected_n_atoms = 9
    expected_n_residues = 1
    expected_n_segments = 1
    expected_n_bonds = 9

    def test_attrs_total_counts(self, top):
        u = mda.Universe(top)
        ag = u.select_atoms("all")
        res = ag.residues
        seg = ag.segments
        assert len(ag) == self.expected_n_atoms
        assert len(res) == self.expected_n_residues
        assert len(seg) == self.expected_n_segments

    def test_aromaticities(self, top, filename):
        expected = np.array([
            atom.IsAromatic() for atom in ob.OBMolAtomIter(filename)])
        assert_equal(expected, top.aromaticities.values)

    def test_elements(self, top, filename):
        expected = np.array([
            OBElementTable().GetSymbol(atom.GetAtomicNum()) for atom in 
            ob.OBMolAtomIter(filename)])
        assert_equal(expected, top.elements.values)

    def test_chiralities(self, top, filename):
        expected = np.array([
            "+" if atom.IsPositiveStereo() else "-" if atom.IsNegativeStereo()
            else "N" for atom in ob.OBMolAtomIter(filename)])
        assert_equal(expected, top.chiralities.values)

    def test_charges(self, top, filename):
        expected = np.array([
            atom.GetPartialCharge() for atom in ob.OBMolAtomIter(filename)])
        assert_allclose(expected, top.charges.values)

    def test_mass_check(self, top, filename):
        expected = np.array([
            atom.GetExactMass() for atom in ob.OBMolAtomIter(filename)])
        assert_allclose(expected, top.masses.values)
