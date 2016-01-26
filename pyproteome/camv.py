"""
This module provides functionality for interacting with CAMV.

Currently limited to importing and outputing scan lists.
"""

# Built-ins
from collections import OrderedDict
import logging
from math import ceil
import os
import subprocess

# Core data analysis libraries
from IPython.display import display
import pandas as pd

from . import utils, modifications


LOGGER = logging.getLogger("pyproteome.camv")
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
CAMV_PATH = utils.which("CAMV.exe")

if CAMV_PATH is None:
    CAMV_PATH = os.path.join(
        THIS_DIR, "..", "..",
        "CAMV", "CAMV", "for_redistribution_files_only", "CAMV.exe",
    )


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
            LOGGER.info(
                "Loading CAMV validation data from \"{}\"".format(
                    os.path.join(*path.split(os.sep)[-2:]),
                )
            )

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
    scan_sets=None,
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
    dict of str, list of int
        Dictionary listing the file names and scans segmented into each file.
    """
    psms = modifications.filter_mod_types(
        psms,
        letter_mod_types=letter_mod_types,
    )

    scan_lists = OrderedDict()

    if len(psms) == 0:
        return scan_lists

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

    if scan_sets is None:
        scan_sets = int(ceil(len(psms) / 2500))

    assert scan_sets >= 0

    slice_sizes = int(ceil(len(scan_list) / scan_sets))

    for i in range(scan_sets):
        out_name = "{}-{}.xls".format(basename, i + 1)
        writer = pd.ExcelWriter(
            os.path.join(
                scan_dir,
                out_name,
            )
        )
        lst = scan_list[i * slice_sizes:(i + 1) * slice_sizes]
        lst = lst.drop_duplicates()
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
    Run the CAMV "Get File" command and "Save Session" commands.

    Parameters
    ----------
    raw_path : str
    xml_path : str
    output_dir : str
    scan_path : str, optional
    save_path : str, optional

    Returns
    -------
    str
    """
    LOGGER.info(
        (
            "Running CAMV import on \"{}\", \"{}\", \"{}\","
            " scans=\"{}\", save=\"{}\""
        ).format(
            os.path.basename(raw_path),
            os.path.basename(xml_path),
            os.path.split(output_dir)[0],
            os.path.basename(scan_path) if scan_path else "",
            os.path.basename(save_path) if save_path else "",
        )
    )
    cmd = [
        CAMV_PATH,
        "import", "true",
        "raw_path", raw_path,
        "xml_path", xml_path,
        "out_path", output_dir,
    ]

    if scan_path:
        cmd += ["sl_path", scan_path]

    if save_path:
        cmd += ["save", "true", "session_path", save_path]

    cmd += ["exit", "true"]

    output = subprocess.check_output(cmd)

    return output


def _run_camv_export(save_path):
    """
    Run CAMV "Load Session" and "Export" commands.

    Parameters
    ----------
    save_path : str

    Returns
    -------
    str
    """
    LOGGER.info(
        "Running CAMV export on \"{}\""
        .format(
            save_path
        )
    )

    cmd = [
        CAMV_PATH,
        "load", "true",
        "session_path", save_path,
        "export", "true",
        "exit", "true",
    ]

    output = subprocess.check_output(cmd)

    return output


def run_camv_validation(scan_lists, force=False):
    """
    Run CAMV on a list of scans.

    Does the initial step of converting, importing, and processing scans.

    Parameters
    ----------
    scan_lists : dict of str, list of int
    force : bool, optional
    """
    for scan_path, scan_list in scan_lists.items():
        # Build a list of paths
        file_name = os.path.basename(scan_path)
        base_name = file_name.rsplit("-", 1)[0]
        base_dir = os.path.abspath("..")
        raw_path = os.path.join(
            base_dir, "MS RAW", base_name + ".raw",
        )
        mascot_xml_path = os.path.join(
            base_dir, "Mascot XMLs", base_name + ".xml",
        )
        camv_output_dir = os.path.join(
            base_dir, "CAMV Output",
        )
        save_path = os.path.join(
            base_dir, "CAMV Sessions", os.path.splitext(file_name)[0] + ".mat"
        )

        if not force and \
           os.path.exists(save_path) and \
           os.stat(save_path).st_size >= 2 ** 12:
            continue

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


def run_camv_export(scan_lists):
    """
    Run CAMV export command.

    Creates a list of excel files for accept / maybe / reject peptides using
    any saved CAMV sessions.

    Parameters
    ----------
    scan_lists : dict of str, list of int
    """
    for scan_path, scan_list in scan_lists.items():
        file_name = os.path.basename(scan_path)
        save_path = os.path.join(
            "..", "CAMV Sessions", os.path.splitext(file_name)[0] + ".mat"
        )
        display(_run_camv_export(save_path))
