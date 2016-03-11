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

from . import paths


PROTEOWIZARD_PATH = None
PROTEOWIZARD_PATH = r"C:\Users\Nader\Dropbox (MIT)\White Lab\CAMV\ProteoWizard\ProteoWizard 3.0.9205"


def fetch_proteowizard():
    pass


def get_scan_data(basename, scans_used):
    raw_path = os.path.join(paths.MS_RAW_DIR, "{}.raw".format(basename))
    ms_convert_path = os.path.join(PROTEOWIZARD_PATH, "msconvert.exe")

    # Create a config file,
    config = tempfile.NamedTemporaryFile(mode="w+", suffix=".txt")
    config.write(
        "filer=\"scanNumber {}\"".format(
            " ".join(str(scan) for scan in scans_used)
        )
    )
    config.flush()

    out_dir = tempfile.mkdtemp()

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
    ms2_data = pymzml.run.Reader(
        out_path,
    )

    prefix = {"mzml": "http://psi.hupo.org/ms/mzml"}
    isolation_mzs = OrderedDict()

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
        print(scan, precursor_scan)
        break

    # del ms2_data
    # shutil.rmtree(out_dir)

    return isolation_mzs, spectrum
