# Testing OpenBabel and Pybel

import MDAnalysis as mda
import openbabel as OB
from openbabel import OBMol, OBConversion
import pybel
from pybel import readfile

import mda_openbabel_converter
from mda_openbabel_converter import OpenBabelParser as OBParser
import pytest  # version 8.2.2
import sys
import numpy as np
from numpy.testing import assert_equal, assert_allclose
# from mda_openbabel_converter.tests.test_data import # something...
from MDAnalysisTests.topology.base import ParserBase  # version 2.7.0

# *** can run with "python -m pytest" but not "pytest" (can't find 
# MDAnalysis) - need to fix this! ***


class OpenBabelParserBase(ParserBase):
    parser = OBParser
    expected_attrs = ['ids', 'names', 'elements', 'masses', 'aromaticities',
                      'resids', 'resnums', 'chiralities',
                      'segids', 'bonds',
                     ]
    
    expected_n_atoms = 0
    expected_n_residues = 1
    expected_n_segments = 1
    expected_n_bonds = 0

    def test_creates_universe(self, filename):
        u = mda.Universe(filename, format='OPENBABEL')
        assert isinstance(u, mda.Universe)

    def test_bonds_total_counts(self, top):
        assert len(top.bonds.values) == self.expected_n_bonds

class TestOpenBabelParserEmpty(OpenBabelParserBase):

    @pytest.fixture(scope='class')
    def filename(self):
        return OBMol()

    expected_attrs = ['ids', 'names', 'elements', 'masses', 'aromaticities',
                      'resids', 'resnums', 'chiralities',
                      'segids', 'bonds',
                     ]
    
    expected_n_atoms = 0
    expected_n_residues = 1
    expected_n_segments = 1
    expected_n_bonds = 0

    @pytest.fixture(scope='class')
    def topology(self, filename):
        topology = OBParser(filename)
        return topology

    def test_attributes(self, topology):

        assert len(np.unique(topology.results.hbonds[:, 0])) == 10
        assert len(topology.results.hbonds) == 32

        reference = {
            'distance': {'mean': 2.7627309, 'std': 0.0905052},
            'angle': {'mean': 158.9038039, 'std': 12.0362826},
        }

        assert_allclose(np.mean(h.results.hbonds[:, 4]),
                        reference['distance']['mean'])
        assert_allclose(np.std(h.results.hbonds[:, 4]),
                        reference['distance']['std'])
        assert_allclose(np.mean(h.results.hbonds[:, 5]),
                        reference['angle']['mean'])
        assert_allclose(np.std(h.results.hbonds[:, 5]),
                        reference['angle']['std'])
        
    def test_atoms_total_counts(self, topology):
        assert len(topology.select_atoms("all")) == self.expected_n_atoms

    def test_residues_total_counts(self, topology):
        assert len(topology.select_atoms("all")) == self.expected_n_atoms

    def test_segments_total_counts(self, topology):
        assert len(topology.select_atoms("all")) == self.expected_n_atoms
        

class TestOpenBabelParserAtomBuild(OpenBabelParserBase):

    @pytest.fixture(scope='class')
    def filename(self):
        return OBMol()

    expected_attrs = ['ids', 'names', 'elements', 'masses', 'aromaticities',
                      'resids', 'resnums', 'chiralities',
                      'segids', 'bonds',
                     ]
    
    # expected_attrs = OpenBabelParserBase.expected_attrs + ['charges']
    
    expected_n_atoms = 2
    expected_n_residues = 1
    expected_n_segments = 1
    expected_n_bonds = 1

    @pytest.fixture(scope='class')
    def parser(self, mol):
        topology = OBParser(mol, **self.kwargs)
        return topology

    def test_hbond_analysis(self, h):

        assert len(np.unique(h.results.hbonds[:, 0])) == 10
        assert len(h.results.hbonds) == 32

        reference = {
            'distance': {'mean': 2.7627309, 'std': 0.0905052},
            'angle': {'mean': 158.9038039, 'std': 12.0362826},
        }

        assert_allclose(np.mean(h.results.hbonds[:, 4]),
                        reference['distance']['mean'])
        assert_allclose(np.std(h.results.hbonds[:, 4]),
                        reference['distance']['std'])
        assert_allclose(np.mean(h.results.hbonds[:, 5]),
                        reference['angle']['mean'])
        assert_allclose(np.std(h.results.hbonds[:, 5]),
                        reference['angle']['std'])
        
mol = OBMol()
print(mol.NumAtoms()) #Should print 0 (atoms)

a = mol.NewAtom()
a.SetAtomicNum(6)   # carbon atom
a.SetVector(0.0, 1.0, 2.0) # coordinates
b = mol.NewAtom()
mol.AddBond(1, 2, 1)   # atoms indexed from 1
print(mol.NumAtoms()) # Should print 2 (atoms)
print(mol.NumBonds()) # Should print 1 (bond)

for i in range(1, mol.NumAtoms()+1):
    atom = mol.GetAtomById(i-1)
    print(a)

# --------

obConversion = OBConversion()
obConversion.SetInAndOutFormats("smi", "mdl")

mol = OBMol()
obConversion.ReadString(mol, "C1=CC=CS1")

#readfile(format="smi", filename="")

print(mol.NumAtoms()) # Should print 5 (atoms)

mol.AddHydrogens()
print(mol.NumAtoms()) # Should print 9 (atoms) after adding hydrogens

outMDL = obConversion.WriteString(mol)
print(outMDL)

# --------


