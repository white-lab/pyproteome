"""
This module provides functionality for interacting with CAMV.

Currently limited to importing and outputing scan lists.
"""

# Built-ins
from collections import OrderedDict
import logging
import os
import subprocess

# Core data analysis libraries
import pandas as pd

from . import utils, modifications


LOGGER = logging.getLogger("pyproteome.camv")
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
CAMV_PATH = os.path.join(THIS_DIR, "CAMV", "CAMV.exe")


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

    scan_list = psms[["First Scan"]]
    scan_list.sort(
        columns=["First Scan"],
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
        scan_lists[out_name] = lst["First Scan"].tolist()
        lst.to_excel(
            writer,
            index=False,
            header=False,
            sheet_name="Scan List",
        )
        writer.save()

    return scan_lists


def _run_camv_get_file(
    raw_path, xml_path, output_dir,
    scan_path=None, save_path=None,
):
    """
    Run the CAMV "Get File" command and then "Save Session".

    Parameters
    ----------
    raw_path : str
    xml_path : str
    output_dir : str
    scan_path : str, optional
    save_path : str, optional
    """
    LOGGER.info(
        "Running CAMV on \"{}\", \"{}\", \"{}\", scans=\"{}\", save=\"{}\""
        .format(
            os.path.basename(raw_path),
            os.path.basename(xml_path),
            os.path.split(output_dir)[0],
            os.path.basename(scan_path) if scan_path else "",
            os.path.basename(save_path) if save_path else "",
        )
    )
    cmd = [
        CAMV_PATH,
        "--get-file",
        raw_path, xml_path, output_dir,
    ]

    if scan_path:
        cmd += scan_path

    if save_path:
        cmd += ["--save-session", save_path]

    output = subprocess.check_call(cmd)

    return output


def run_camv_validation(scan_lists):
    """
    Run CAMV on a list of scans.

    Does the initial step of converting, importing, and processing scans.

    Parameters
    ----------
    scan_lists : dict of str, list of str
    """
    for scan_path, scan_list in scan_lists.items():
        # Build a list of paths
        file_name = os.path.basename(scan_path)
        base_name = file_name.rsplit("-", 1)[0]
        raw_path = os.path.join(
            "..", "MS RAW", base_name + ".raw",
        )
        mascot_xml_path = os.path.join(
            "..", "Mascot XMLs", base_name + ".xml",
        )
        camv_output_dir = os.path.join(
            "..", "CAMV Output",
        )
        save_path = os.path.join(
            "..", "CAMV Sessions", os.path.splitext(file_name)[0] + ".mat"
        )

        # Run CAMV
        _run_camv_get_file(
            raw_path,
            mascot_xml_path,
            camv_output_dir,
            scan_path=scan_path,
            save_path=save_path,
        )

        # Check save files exists and is >= 4 MB
        if not os.path.exists(save_path) or \
           os.stat(save_path).st_size < 2 ** 12:
            raise Exception(
                "CAMV did not create a save file for {}".format(save_path)
            )
