"""
pyproteome is a Python package for interacting with proteomics data.

It includes modules for loading, processing, and analyzing data by mass
spectometry. Currently it only supports data produced by MASCOT / Discoverer
and CAMV.
"""

from .version import __version__

from . import (
    analysis, bca, data_sets, discoverer, icelogo, levels, loading, logo,
    modification, motif, paths, plogo, pride, protein, sequence, utils,
    version, volcano, weblogo,
)

__all__ = [
    "analysis",
    "bca",
    "data_sets",
    "discoverer",
    "icelogo",
    "levels",
    "loading",
    "logo",
    "modification",
    "motif",
    "paths",
    "plogo",
    "pride",
    "protein",
    "sequence",
    "utils",
    "version",
    "volcano",
    "weblogo",
    "__version__",
]
