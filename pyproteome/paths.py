"""
This module tracks the path to user data files. Developers can override paths
here when using a custom data hierarchy.
"""

import logging
import os

LOGGER = logging.getLogger("pyproteome.paths")

BASE_DIR = os.path.abspath("..")
BASE_DIR_OPTS = (
    os.path.abspath("."),
    os.path.abspath(".."),
)

(
    BCA_ASSAY_DIR,
    CAMV_OUT_DIR,
    CAMV_SESS_DIR,
    MS_SEARCHED_DIR,
    MS_RAW_DIR,
    SCRIPTS_DIR,
    FIGURES_DIR,
) = (None,) * 7

BCA_NAME = "BCA Protein Assays"
CAMV_NAME = "CAMV Output"
CAMV_SESS_NAME = "CAMV Sessions"
MS_SEARCHED_NAME = "Searched"
MS_RAW_NAME = "MS RAW"
SCRIPTS_NAME = "Script"
FIGURES_NAME = "Figures"

DIR_NAMES = (
    BCA_NAME,
    CAMV_NAME,
    CAMV_SESS_NAME,
    MS_SEARCHED_NAME,
    MS_RAW_NAME,
    SCRIPTS_NAME,
    FIGURES_NAME,
)


def set_base_dir(path):
    """
    Set the base directory containing the search / raw / scripts / figures
    folders.

    Parameters
    ----------
    path : str
    """
    global \
        BCA_ASSAY_DIR, CAMV_OUT_DIR, CAMV_SESS_DIR, MASCOT_XML_DIR, \
        MS_SEARCHED_DIR, MS_RAW_DIR, SCRIPTS_DIR, FIGURES_DIR

    BCA_ASSAY_DIR = os.path.join(path, BCA_NAME)
    CAMV_OUT_DIR = os.path.join(path, CAMV_NAME)
    CAMV_SESS_DIR = os.path.join(path, CAMV_SESS_NAME)
    MS_SEARCHED_DIR = os.path.join(path, MS_SEARCHED_NAME)
    MS_RAW_DIR = os.path.join(path, MS_RAW_NAME)
    SCRIPTS_DIR = os.path.join(path, SCRIPTS_NAME)
    FIGURES_DIR = os.path.join(path, FIGURES_NAME)


def find_base_dir():
    """
    Finds the base directory containing the search / raw / scripts / figures
    folders. May be the current working directory or a parent of it.

    Returns
    -------
    path : str
    """
    for opt in BASE_DIR_OPTS:
        if any(
            os.path.exists(os.path.join(opt, i))
            for i in DIR_NAMES
        ):
            return opt

    LOGGER.warning(
        "Unable to find any expected data directories relative to {}: {}"
        .format(
            os.getcwd(),
            ", ".join(DIR_NAMES),
        )
    )
    LOGGER.warning(
        "Setting base path to {}, consider calling paths.set_base_dir()"
        .format(os.getcwd())
    )

    return BASE_DIR_OPTS[0]


set_base_dir(find_base_dir())
