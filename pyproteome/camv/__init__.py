"""
This module provides functionality for interacting with CAMV.

Currently limited to importing and outputing scan lists.
"""

from __future__ import absolute_import, division

import logging
import os

import pandas as pd

import pyproteome as pyp


try:
    FileNotFoundError
except NameError:
    FileNotFoundError = (IOError, OSError)

LOGGER = logging.getLogger("pyproteome.camv")
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
CAMV_PATH = pyp.utils.which("CAMV.exe")

if CAMV_PATH is None:
    CAMV_PATH = os.path.abspath(
        os.path.join(
            THIS_DIR, "..",
            "CAMV", "CAMV", "for_redistribution_files_only", "CAMV.exe",
        )
    )


def load_camv_validation(basename):
    """
    Load validation data produced by CAMV.

    Parameters
    ----------
    basename : str

    Returns
    -------
    accepted : :class:`pandas.DataFrame`
    maybed : :class:`pandas.DataFrame`
    rejected : :class:`pandas.DataFrame`
    """
    accepted = None
    maybed = None
    rejected = None

    def _try_open_xls(path, existing=None):
        try:
            df = pd.read_csv(path, sep="\t")
        except OSError:
            return existing
        else:
            LOGGER.info(
                "Loading CAMV validation data from \"{}\"".format(
                    os.path.join(*path.split(os.sep)[-2:]),
                )
            )
            in_name = os.path.splitext(os.path.basename(path))[0]
            df["Scan Paths"] = pd.Series(
                [set([in_name])] * len(df.index)
            )

            if existing is not None:
                df = pd.concat([existing, df])

            return df

    try:
        files = os.listdir(pyp.paths.CAMV_OUT_DIR)
    except FileNotFoundError:
        return accepted, maybed, rejected

    for filename in files:
        if filename.startswith(basename):
            base_dir = os.path.join(pyp.paths.CAMV_OUT_DIR, filename)

            accept_path = os.path.join(base_dir, "accept.xls")
            maybe_path = os.path.join(base_dir, "maybe.xls")
            reject_path = os.path.join(base_dir, "reject.xls")

            accepted = _try_open_xls(accept_path, existing=accepted)
            maybed = _try_open_xls(maybe_path, existing=maybed)
            rejected = _try_open_xls(reject_path, existing=rejected)

    return accepted, maybed, rejected
