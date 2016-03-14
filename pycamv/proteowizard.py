"""
Provides functions for interacting with MS data through ProteoWizard.
"""

from collections import OrderedDict
import os
import re
import shutil
import subprocess
import tempfile

import pymzml

from pyproteome import paths


PROTEOWIZARD_PATH = None
PROTEOWIZARD_PATH = (
    r"C:\Users\Nader\Dropbox (MIT)\White Lab\CAMV\ProteoWizard"
    r"\ProteoWizard 3.0.9205"
)


def fetch_proteowizard():
    pass


def _raw_to_mzml(basename, scans, out_dir, mz_window=None):
    """
    Parameters
    ----------
    basename : str
    scans : list of int
    out_dir : str
    mz_window : list of int, optional
    """
    ms_convert_path = os.path.join(PROTEOWIZARD_PATH, "msconvert.exe")
    raw_path = os.path.join(paths.MS_RAW_DIR, "{}.raw".format(basename))

    # Create a config file,
    config = tempfile.NamedTemporaryFile(mode="w+", suffix=".txt")

    config.write(
        "filter=\"scanNumber {}\"\n".format(
            " ".join(str(scan) for scan in scans)
        )
    )

    if mz_window:
        config.write(
            "filter=\"mzWindow [{},{}]\"\n".format(
                mz_window[0], mz_window[1],
            )
        )

    config.flush()

    cmd = [
        ms_convert_path,
        raw_path,
        "-o", out_dir,
        "--mzML",
        "-c", config.name,
    ]

    subprocess.check_call(cmd)
    config.close()

    out_path = os.path.join(out_dir, "{}.mzML".format(basename))
    data = pymzml.run.Reader(
        out_path,
    )

    return data


def get_scan_data(basename, scans_used):
    out_dir = tempfile.mkdtemp()
    ms2_data = _raw_to_mzml(basename, scans_used, out_dir)

    prefix = {"mzml": "http://psi.hupo.org/ms/mzml"}
    isolation_mzs = OrderedDict()
    precursor_scans = set()

    for spectrum in ms2_data:
        if spectrum["ms level"] != 2:
            continue

        scan = spectrum["id"]
        isolation_mz = spectrum["selected ion m/z"]
        isolation_mzs[scan] = isolation_mz

        spectrum_ref = spectrum.xmlTreeIterFree.find(
            "mzml:precursorList/mzml:precursor", prefix,
        ).get("spectrumRef")
        precursor_scan = re.search("scan=(\d+)", spectrum_ref).group(1)
        precursor_scans.add(precursor_scan)

    ms_data = _raw_to_mzml(basename, sorted(precursor_scans), out_dir)

    # del ms2_data
    # shutil.rmtree(out_dir)

    return isolation_mzs, spectrum
