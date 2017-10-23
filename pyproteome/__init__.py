"""
pyproteome is a Python package for interacting with proteomics data.

It includes modules for loading, processing, and analyzing data by mass
spectometry. Currently it only supports data produced by MASCOT / Discoverer
and CAMV.
"""

from .version import __version__

from . import (
    analysis, bca, data_sets, discoverer, levels, loading,
    modification, paths, pride, protein, sequence, utils,
    version, volcano,
)
from .motifs import (
    logo, icelogo, motif, neighborhood, plogo, phosphosite, weblogo,
)

__all__ = [
    "analysis",
    "bca",
    "data_sets",
    "discoverer",
    "levels",
    "loading",
    "modification",
    "paths",
    "pride",
    "protein",
    "sequence",
    "utils",
    "version",
    "volcano",
    "__version__",

    "icelogo",
    "logo",
    "motif",
    "neighborhood",
    "plogo",
    "phosphosite",
    "weblogo",
]
