'''
This module tracks the path to user data files. Developers can override paths
here when using a custom data hierarchy.
'''

import logging
import os

LOGGER = logging.getLogger('pyproteome.paths')

BASE_DIR = os.path.abspath('.')
'''
Location of the base directory containing proteomics data.
By default this is set to the current or parent directory, whichever contains
any folders matching the expected directory structure.
'''

BASE_DIR_OPTS = (
    os.path.abspath('.'),
    os.path.abspath('..'),
)

CAMV_OUT_DIR = None
'''
Location of the directory containing validated CAMV data.
By default it is set to :const:`.FIGURES_NAME` in the current or parent directory.
'''
MS_SEARCHED_DIR = None
'''
Location of the directory containing Proteome Discoverer .msf search files.
By default it is set to :const:`.FIGURES_NAME` in the current or parent directory.
'''
MS_RAW_DIR = None
'''
Location of the directory containing raw mass spectrometry files.
By default it is set to :const:`.FIGURES_NAME` in the current or parent directory.
'''
FIGURES_DIR = None
'''
Location of the directory for saving output figures.
By default it is set to :const:`.FIGURES_NAME` in the current or parent directory.
'''

CAMV_NAME = 'CAMV Output'
'''
Name of the directory containing validated CAMV data.
'''

MS_SEARCHED_NAME = 'Searched'
'''
Name of the directory containing Proteome Discoverer .msf search files.
'''

MS_RAW_NAME = 'MS RAW'
'''
Name of the directory containing raw mass spectrometry files.
'''

FIGURES_NAME = 'Figures'
'''
Name of the directory for saving output figures.
'''

DIR_NAMES = (
    CAMV_NAME,
    MS_SEARCHED_NAME,
    MS_RAW_NAME,
    FIGURES_NAME,
)


def set_base_dir(path):
    '''
    Set the base directory containing the search / raw / figures
    folders.

    Parameters
    ----------
    path : str
    '''
    global \
        CAMV_OUT_DIR, MS_SEARCHED_DIR, MS_RAW_DIR, \
        FIGURES_DIR

    CAMV_OUT_DIR = os.path.join(path, CAMV_NAME)
    MS_SEARCHED_DIR = os.path.join(path, MS_SEARCHED_NAME)
    MS_RAW_DIR = os.path.join(path, MS_RAW_NAME)
    FIGURES_DIR = os.path.join(path, FIGURES_NAME)


def find_base_dir():
    '''
    Finds the base directory containing the search / raw / scripts / figures
    folders. May be the current working directory or a parent of it.

    Returns
    -------
    path : str
    '''
    for opt in BASE_DIR_OPTS:
        if any(
            os.path.exists(os.path.join(opt, i))
            for i in DIR_NAMES
        ):
            return opt

    LOGGER.warning(
        'Unable to find any expected data directories relative to {}: {}'
        .format(
            os.getcwd(),
            ', '.join(DIR_NAMES),
        )
    )
    LOGGER.warning(
        'Setting base path to {}, consider calling paths.set_base_dir()'
        .format(os.getcwd())
    )

    return BASE_DIR_OPTS[0]


set_base_dir(find_base_dir())
