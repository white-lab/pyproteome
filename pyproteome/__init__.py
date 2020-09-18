
from .utils import DEFAULT_DPI
from .analysis import (
    correlation, tables, volcano,
)
from .motifs import (
    logo, motif, phosphosite,
)
from . import (
    analysis, data_sets, discoverer, levels, loading, motifs, paths,
    pride, species, utils, version,
)
from . import cluster, pathways

try:
    from IPython import get_ipython
    from IPython.core.magic import register_line_magic
except ImportError:
    get_ipython = None
    register_line_magic = None


def import_all(line=None):
    '''
    Inialize and import many packages using IPython Notebooks magic.

    Imports numpy pandas, seaborn sklearn, and pyproteome packages.
    Sets visual display options for matplotlib and adds a logging handlers.
    Also applies auto-reload to pyproteome for developers.

    Examples
    --------
    >>> from pyproteome import *
    >>> %import_all
    '''
    if get_ipython is None:
        return

    ip = get_ipython()
    ip.run_line_magic(
        'config',
        'InlineBackend.figure_formats = [\'retina\']',
    )
    ip.run_line_magic('load_ext', 'autoreload')
    ip.run_line_magic('autoreload', '2')
    ip.run_line_magic('aimport', 'pyproteome')
    ip.run_line_magic('aimport', 'brainrnaseq')
    ip.run_line_magic('pylab', 'inline')

    ip.ex(
        '\n'.join([
            'from collections import OrderedDict, Counter',
            'import logging',
            'import os',
            'import pickle',
            'from IPython.display import display, SVG, Image',

            'import numpy as np',
            'import pandas as pd',
            'import seaborn as sns',
            'import sklearn',
            'import pyproteome as pyp',
            'import brainrnaseq as brs',

            'pylab.rcParams[\'figure.figsize\'] = (6, 4)',
            'pylab.rcParams[\'figure.dpi\'] = 144',
            'pylab.rcParams[\'mathtext.default\'] = \'regular\'',
            'pylab.rcParams[\'figure.max_open_warning\'] = 0',

            'sns.set_style(\'white\')',
            'sns.set_context(\'notebook\')',

            'pd.set_option(\'display.max_colwidth\', 500)',
            'pd.set_option(\'display.max_rows\', 500)',

            'formatter = logging.Formatter(\'%(asctime)s\t%(name)s\t'
            '%(levelname)s\t%(message)s\', datefmt=\'%I:%M:%S %p\')',

            'root = logging.getLogger()',

            'if not root.handlers: '
            'handler = logging.StreamHandler(); '
            'handler.setFormatter(formatter); '
            'root.setLevel(logging.INFO); '
            'root.addHandler(handler)',
        ])
    )


if (
    register_line_magic is not None and
    get_ipython is not None and
    get_ipython() is not None
):
    import_all = register_line_magic(import_all)

__all__ = [
    'analysis',
    'cluster',
    'correlation',
    'data_sets',
    'discoverer',
    'levels',
    'loading',
    'motifs',
    'phosphosite',
    'paths',
    'pathways',
    'pride',
    'species',
    'tables',
    'utils',
    'version',
    'volcano',

    'import_all',
    'DEFAULT_DPI',

    'logo',
    'motif',
]
