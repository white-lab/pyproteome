"""
This module provides functionality for accessing public proteome data through
PRIDE / Proteome Xchange.
"""

from io import BytesIO
import os
# XXX: This should be a safer alternative package. Otherwise users could be
# DOS-d by a MITM attack
import xml.etree.ElementTree as ET
from zipfile import ZipFile

import requests

META_DATA_URL = (
    "http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID={}&"
    "outputMode=XML&test=no"
)


def list_data_set(accession):
    """
    Lists files contained in a deposition on PRIDE.

    Information is fetched from pride.META_DATA_URL.

    Parameters
    ----------
    accession : str

    Returns
    -------
    list of xml.etree.ElementTree
    """
    assert accession.startswith("PXD")

    proj_id = int(accession[3:])

    # Fetch xml file using requests
    meta_data = requests.get(META_DATA_URL.format(proj_id))
    meta_data.raise_for_status()

    root = ET.fromstring(meta_data.text)

    return root.findall("DatasetFileList/DatasetFile")


def fetch_data_set(accession, unzip=None):
    """
    Fetches files from a deposition on PRIDE.

    Parameters
    ----------
    accession : str
        A PRIDE accession ID. i.e. "PXD001038"
    unzip : dict of str, str, optional
        Maps zip files to the directories they should be extracted into. In
        unset, all zip files for a project will be downloaded and extracted
        into the current working directory.

    Examples
    --------
    >>> from pyproteome import pride
    >>> pride.fetch_data_set(
    ...     "PXD001038",
    ...     unzip={"HJ070512_OCTFF_B2_All5Fractions_PeptideSummary.zip": "."},
    ... )
    """
    files = list_data_set(accession)

    for file_root in files:
        if unzip and file_root.get("name") not in unzip:
            continue

        folder = unzip[file_root.get("name")] if unzip else os.getcwd()

        file_url = file_root.find("cvParam").get("value")

        if not file_url.endswith(".zip"):
            raise Exception("Unknown file type: {}".format(file_url))

        # Requests cannot fetch FTP files
        if file_url.startswith("ftp://"):
            file_url = "http://{}".format(file_url[len("ftp://"):])

        zip_file = requests.get(file_url, stream=True)
        zip_contents = ZipFile(BytesIO(zip_file.content))

        for name in zip_contents.namelist():
            zip_contents.extract(name, folder)
