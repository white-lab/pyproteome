"""
pyproteome is a Python package for interacting with proteomics data.

It includes modules for loading, processing, and analyzing data by mass
spectometry. Currently it only supports data produced by MASCOT / Discoverer
and CAMV.
"""

from . import (
    analysis, bca, cluster, data_sets, discoverer, levels, loading,
    modification, paths, pride, protein, sequence, tables, utils,
    version, volcano,
)
from .motifs import (
    logo, icelogo, motif, neighborhood, plogo, phosphosite, weblogo,
)

from IPython import get_ipython
from IPython.core.magic import register_line_magic


if get_ipython() is not None:
    @register_line_magic
    def import_all(line=None):
        ip = get_ipython()
        ip.run_line_magic(
            "config",
            "InlineBackend.figure_formats = ['retina']",
        )
        ip.run_line_magic("load_ext", "autoreload")
        ip.run_line_magic("autoreload", "2")
        ip.run_line_magic("aimport", "pyproteome")
        ip.run_line_magic("aimport", "brainrnaseq")
        ip.run_line_magic("pylab", "inline")

        ip.ex(
            "\n".join([
                "from collections import OrderedDict, Counter",
                "import os",
                "import pickle",
                "from IPython.display import display",

                "import numpy as np",
                "import pandas as pd",
                "import seaborn as sns",
                "import sklearn",

                "pylab.rcParams['figure.figsize'] = (12, 8)",
                "pylab.rcParams['mathtext.default'] = 'regular'",
                "pylab.rcParams['figure.max_open_warning'] = 0",

                "sns.set_style('white')",
                "sns.set_context('notebook')",

                'pd.set_option("display.max_colwidth", 500)',
                'pd.set_option("display.max_rows", 500)',

                "import logging",
                "root = logging.getLogger()",
                "if not root.handlers: root.addHandler(logging.StreamHandler())",
                "logging.getLogger().setLevel(logging.INFO)",
            ])
        )


__all__ = [
    "analysis",
    "bca",
    "cluster",
    "data_sets",
    "discoverer",
    "levels",
    "loading",
    "modification",
    "paths",
    "pride",
    "protein",
    "sequence",
    "tables",
    "utils",
    "version",
    "volcano",
    "import_all",

    "icelogo",
    "logo",
    "motif",
    "neighborhood",
    "plogo",
    "phosphosite",
    "weblogo",
]
