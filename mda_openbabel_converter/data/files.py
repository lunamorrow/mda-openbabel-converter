"""
Location of data files
======================

Use as ::

    from mda_openbabel_converter.data.files import *

"""

__all__ = [
    "MDANALYSIS_LOGO",  # example file of MDAnalysis logo
]

import importlib.resources
MDANALYSIS_LOGO = importlib.resources.files(__name__) / "mda.txt"
