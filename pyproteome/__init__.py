"""
pyproteome is a Python package for interacting with proteomics data.

It includes modules for loading, processing, and analyzing data by mass
spectometry. Currently it only supports data produced by MASCOT / Discoverer
and CAMV.
"""

from .version import __version__
from .modification import Modification, Modifications
from .sequence import Sequence, ProteinMatch
from .motif import Motif
from .protein import Protein, Proteins
from .data_sets import DataSet
