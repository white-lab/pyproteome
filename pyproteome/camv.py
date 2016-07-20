"""
This module provides functionality for interacting with CAMV.

Currently limited to importing and outputing scan lists.
"""

# Built-ins
from __future__ import division

from collections import OrderedDict
import logging
from math import ceil
import os
import subprocess

# Core data analysis libraries
from IPython.display import display
import pandas as pd

from . import modification, paths, utils


LOGGER = logging.getLogger("pyproteome.camv")
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
CAMV_PATH = utils.which("CAMV.exe")

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

    for filename in os.listdir(paths.CAMV_OUT_DIR):
        if filename.startswith(basename):
            base_dir = os.path.join(paths.CAMV_OUT_DIR, filename)

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
    psms : :class:`pandas.DataFrame`
    basename : str, optional
    letter_mod_types : list of tuple of str or None, str or None, optional
    scan_sets : int, optional

    Returns
    -------
    :class:`pandas.DataFrame`
    dict of str, list of int
        Dictionary listing the file names and scans segmented into each file.
    """
    filtered = modification.filter_mod_types(
        psms,
        letter_mod_types=letter_mod_types,
    )

    scan_lists = OrderedDict()

    if len(filtered) == 0:
        return psms, scan_lists

    psms = psms.copy()
    psms["Scan Paths"] = pd.Series([set("")] * len(psms.index))

    scan_list = filtered.sort(
        columns=["Protein Group Accessions", "First Scan"],
        ascending=True,
    )[["First Scan"]]

    scan_list.drop_duplicates(
        subset="First Scan",
        inplace=True,
    )

    # Export as XLS
    utils.make_folder(paths.SCAN_LISTS_DIR)

    if scan_sets is None:
        scan_sets = int(ceil(len(scan_list) / 1000))
        LOGGER.info(
            "Splitting {} unique scans into {} lists".format(
                len(scan_list),
                scan_sets,
            )
        )

    assert scan_sets >= 0

    slice_sizes = int(ceil(len(scan_list) / scan_sets))

    for i in range(scan_sets):
        out_name = "{}-{}.xls".format(basename, i + 1)
        writer = pd.ExcelWriter(
            os.path.join(
                paths.SCAN_LISTS_DIR,
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

        # XXX: Better way to do this?
        for index, row in psms[
            psms["First Scan"].isin(scan_lists[out_name])
        ].iterrows():
            psms.set_value(
                index,
                "Scan Paths",
                set([os.path.splitext(out_name)[0]]),
            )

    return psms, scan_lists


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
        cmd += [
            "save", "true",
            "session_path", save_path,
        ]

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


def _run_each_camv_validation(scan_path, scan_list, force=False):
    # Build a list of paths
    scan_path = os.path.join(
        paths.SCAN_LISTS_DIR, scan_path,
    )
    file_name = os.path.basename(scan_path)
    base_name = file_name.rsplit("-", 1)[0]
    raw_path = os.path.join(
        paths.MS_RAW_DIR, base_name + ".raw",
    )
    mascot_xml_path = os.path.join(
        paths.MASCOT_XML_DIR, base_name + ".xml",
    )
    save_path = os.path.join(
        paths.CAMV_SESS_DIR, os.path.splitext(file_name)[0] + ".mat"
    )

    # Skip archived scan lists
    for store_ext in [".7z", ".zip", ".tar", ".tar.gz", ".tar.bz2"]:
        if os.path.exists(save_path + store_ext) or \
           os.path.exists(os.path.splitext(save_path)[0] + store_ext):
            return

    if not force and \
       os.path.exists(save_path) and \
       os.stat(save_path).st_size >= 2 ** 12:
        return

    # Run CAMV
    _run_camv_get_file(
        raw_path,
        mascot_xml_path,
        paths.CAMV_OUT_DIR,
        scan_path=scan_path,
        save_path=save_path,
    )

    # Check save files exists and is >= 4 MB
    if not os.path.exists(save_path) or \
       os.stat(save_path).st_size < 2 ** 12:
        raise Exception(
            "CAMV did not create a save file for {}".format(save_path)
        )


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
        _run_each_camv_validation(scan_path, scan_list, force=force)


def run_camv_export(scan_lists=None):
    """
    Run CAMV export command.

    Creates a list of excel files for accept / maybe / reject peptides using
    any saved CAMV sessions. If scan_lists is None, everything in the "CAMV
    Sessions" directory will be exported.

    Parameters
    ----------
    scan_lists : dict of str, list of int, optional
    """
    if scan_lists is None:
        scan_paths = [
            path
            for path in os.listdir(paths.CAMV_SESS_DIR)
            if path.endswith(".mat")
        ]
    else:
        scan_paths = list(scan_lists.keys())

    for scan_path in scan_paths:
        file_name = os.path.basename(scan_path)
        save_path = os.path.join(
            paths.CAMV_SESS_DIR, os.path.splitext(file_name)[0] + ".mat"
        )
        display(_run_camv_export(save_path))
