"""
Location of data files
======================

Use as ::

    from mda_openbabel_converter.data.files import *

"""

__all__ = [
    "MDANALYSIS_LOGO",  # example file of MDAnalysis logo
]

from importlib.resources import files
MDANALYSIS_LOGO = files("mda_openbabel_converter") / "data" / "mda.txt"
CRN = files("mda_openbabel_converter") / "data" / "1crn.pdb"