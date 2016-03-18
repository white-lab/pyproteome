"""
Provides functions for interacting with MS data through ProteoWizard.
"""

import logging
import os
import platform
import requests
import shutil
import subprocess
import sys
import tempfile

import pymzml

from pyproteome import paths


LOGGER = logging.getLogger("pycamv.proteowizard")


THIS_DIR = os.path.abspath(os.path.dirname(__file__))
PROTEOWIZARD_DIR = os.path.join(THIS_DIR, "ProteoWizard")

PROTEOWIZARD_VERSION = "3.0.9393"
PROTEOWIZARD_PATH = os.path.join(
    PROTEOWIZARD_DIR,
    "ProteoWizard {}".format(PROTEOWIZARD_VERSION),
)

PROTEOWIZARD_MSI_32_URL = (
    "https://www.dropbox.com/s/5q7flvb32vgdh7d/"
    "pwiz-setup-3.0.9393-x86.msi?dl=1"
)

PROTEOWIZARD_MSI_64_URL = (
    "https://www.dropbox.com/s/7bz3wpgjj7zd99a/"
    "pwiz-setup-3.0.9393-x86_64.msi?dl=1"
)


def fetch_proteowizard(url=None):
    if os.path.exists(PROTEOWIZARD_PATH):
        return

    LOGGER.info("ProteoWizard not installed, fetching now.")

    if platform.system() not in ["Windows"]:
        raise Exception("Proteowizard install not supported on your platform")

    if url is None:
        if platform.architecture()[0] == "64bit":
            url = PROTEOWIZARD_MSI_64_URL
        else:
            url = PROTEOWIZARD_MSI_32_URL

    tmpdir = tempfile.mkdtemp()
    out_path = os.path.join(tmpdir, url.rsplit("/", 1)[1].rsplit("?")[0])

    # Download the .msi file
    with open(out_path, mode="wb") as f:
        response = requests.get(url, stream=True)

        if not response.ok:
            raise Exception("Unable to download file: {}".format(response))

        for block in response.iter_content(1024):
            f.write(block)

    # Extract the msi file's contents
    extract_path = os.path.join(tmpdir, "msi_extract")
    cmd = [
        "msiexec",
        "/a",
        out_path,
        "/qb",
        "TARGETDIR=\"{}\"".format(extract_path),
    ]
    subprocess.check_call(" ".join(cmd), shell=True)

    # Copy the msi file's contents to PROTEOWIZARD_DIR
    src = os.path.join(
        extract_path,
        "PFiles",
        "ProteoWizard",
    )

    shutil.rmtree(PROTEOWIZARD_DIR)
    shutil.copytree(src, PROTEOWIZARD_DIR)
    shutil.rmtree(tmpdir)


def raw_to_mzml(basename, out_dir, scans=None, mz_window=None):
    """
    Covert a RAW file to .mzML using ProteoWizard.

    Parameters
    ----------
    basename : str
    out_dir : str
    scans : list of int, optional
    mz_window : list of int, optional

    Returns
    -------
    :class:`pymzml.run.Reader<run.Reader>`
    """
    fetch_proteowizard()

    ms_convert_path = os.path.join(PROTEOWIZARD_PATH, "msconvert.exe")
    raw_path = os.path.join(paths.MS_RAW_DIR, "{}.raw".format(basename))

    # Create a config file,
    config, config_path = tempfile.mkstemp(suffix=".txt", text=True)

    with os.fdopen(config, "w+") as config:
        if scans:
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

    # Run msconvert to convert raw file to mzML
    LOGGER.info("Converting \"{}\" to .mzML format.".format(raw_path))

    cmd = [
        ms_convert_path,
        raw_path,
        "-o", out_dir,
        "--mzML",
        "-c", config_path,
    ]

    out = subprocess.check_output(cmd)
    LOGGER.debug(out.decode(sys.stdout.encoding))

    os.remove(config_path)

    # Read the file into memory using pymzml
    out_path = os.path.join(out_dir, "{}.mzML".format(basename))
    data = pymzml.run.Reader(
        out_path,
        extraAccessions=[
            ("MS:1000827", ["value"]),  # isolation window target m/z
            ("MS:1000828", ["value"]),  # isolation window lower offset
            ("MS:1000829", ["value"]),  # isolation window upper offset
            ("MS:1000512", ["value"]),  # filter string
        ],
    )

    return data
