"""
Documentation...
"""

import MDAnalysis as mda
from MDAnalysis.topology.base import TopologyReaderBase

class OpenBabelParser():
    """
    Inherits from TopologyReaderBase and converts an OpenBabel OBMol to a 
    MDAnalysis Topology or adds it to a pre-existing Topology. This parser will 
    does not work in the reverse direction.
    """
    def _format_hint(thing):
        """
        Base func8on to check if the parser can actually parse this “thing” 
        (i.e., is it a valid OpenBabel OBMol with no missing information, that 
        can be converted to a MDAnalysis Topology?)
        """
        pass

    def parse(self, **kwargs):
        """
        Accepts an OpenBabel OBMol and returns a MDAnalysis Topology. Will need
        to extract the number of atoms, number of residues, number of segments,
        atom_residue index, residue_segment index and other attributes from the
        OBMol to initialise a new Topology. 
        """
        pass