"""
This module provides functionality for interacting with CAMV.

Currently limited to importing and outputing scan lists.
"""

# Built-ins
from collections import OrderedDict
import logging
import os

# Core data analysis libraries
import pandas as pd

from . import utils, modifications


LOGGER = logging.getLogger("pyproteome.camv")


def load_camv_validation(basename):
    """
    Load validation data produced by CAMV.

    Parameters
    ----------
    basename : str

    Returns
    -------
    accepted : pandas.DataFrame
    maybed : pandas.DataFrame
    rejected : pandas.DataFrame
    """
    accepted = None
    maybed = None
    rejected = None

    def _try_open_xls(path, existing=None):
        try:
            df = pd.read_csv(path, sep="\t")
        except OSError:
            return None
        else:
            LOGGER.info("Loading CAMV validation data from {}".format(path))

            if existing is not None:
                df = pd.concat([existing, df])

            return df

    camv_dir = os.path.join(
        "..", "CAMV Output",
    )

    for filename in os.listdir(camv_dir):
        if filename.startswith(basename):
            base_dir = os.path.join(camv_dir, filename)

            accept_path = os.path.join(base_dir, "accept.xls")
            maybe_path = os.path.join(base_dir, "maybe.xls")
            reject_path = os.path.join(base_dir, "reject.xls")

            accepted = _try_open_xls(accept_path, existing=accepted)
            maybed = _try_open_xls(maybe_path, existing=maybed)
            rejected = _try_open_xls(reject_path, existing=rejected)

    return accepted, maybed, rejected


def output_scan_list(
    psms,
    basename="Scan List",
    letter_mod_types=None,
    scan_sets=1,
):
    """
    Write a list of scans to file.

    Can select for any peptides a given amino acid / modification type
    (i.e. letter_mod_types=[("Y", "Phospho")] ).

    Parameters
    ----------
    psms : pandas.DataFrame
    basename : str, optional
    letter_mod_types : list of tuple of str or None, str or None, optional
    scan_sets : int, optional

    Returns
    -------
    pandas.DataFrame
        Scan list that is also saved to file
    dict of str, list of str
        Dictionary listing the file names and scans segmented into each file.
    """
    assert scan_sets >= 0

    psms = modifications.filter_mod_types(
        psms,
        letter_mod_types=letter_mod_types,
    )

    scan_list = psms["First Scan"]
    scan_list.sort(
        ascending=True,
        inplace=True,
    )

    # Export as XLS
    scan_dir = os.path.join(
        "..", "Scan Lists",
    )

    utils.make_folder(scan_dir)

    slice_sizes = len(scan_list) // scan_sets + (1)

    scan_lists = OrderedDict()

    for i in range(scan_sets):
        out_name = "{}-{}.xls".format(basename, i + 1)
        writer = pd.ExcelWriter(
            os.path.join(
                scan_dir,
                out_name,
            )
        )
        lst = scan_list[i * slice_sizes:(i + 1) * slice_sizes]
        scan_lists[out_name] = lst.tolist()
        lst.to_excel(
            writer,
            index=False,
            header=False,
            sheet_name="Scan List",
        )
        writer.save()

    return scan_lists
