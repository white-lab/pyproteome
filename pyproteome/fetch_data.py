"""
This module provides functionality for fetching protein data from UniProt.

Caches fetched protein data for faster re-use.
"""

# Built-ins
import os
import re
import shutil
import sqlite3
import tempfile

# Core data analysis libraries
import pandas as pd

import uniprot

from . import paths


RE_ACCESSION = re.compile("\[([A-Za-z0-9]+_[A-Z]+)\]")
RE_DISCOVERER_ACCESSION = re.compile(r"^>sp\|([\dA-Za-z]+)\|[\dA-Za-z_]+ .*$")

UNIPROT_DATA = {}


def fetch_uniprot_data(accessions):
    """
    Fetch UniProt protein descriptions, gene names, sequences, etc.

    All information is stored in UNIPROT_DATA and can be accessed with
    get_uniprot_data().

    Parameters
    ----------
    accessions : list of str
    """
    accessions = set(accessions).difference(UNIPROT_DATA)

    if not accessions:
        return

    cache_dir = tempfile.mkdtemp(suffix="uniprot")

    UNIPROT_DATA.update(
        uniprot.get_metadata_with_some_seqid_conversions(
            accessions,
            cache_dir=cache_dir,
        )
    )

    shutil.rmtree(cache_dir)


def prefetch_all_uniprot():
    """
    Fetch data for all accesions found in MS Searched directory.

    Pulls all UniProt accession IDs from all "_psms.txt" files.
    """
    accessions = set()

    for filename in os.listdir(paths.MS_SEARCHED_DIR):
        if not filename.endswith("_psms.txt"):
            continue

        psms = pd.read_table(os.path.join(paths.MS_SEARCHED_DIR, filename))
        psms.dropna(
            subset=["Protein Group Accessions"],
            inplace=True,
        )
        accessions.update(
            acc.strip()
            for row_str in psms["Protein Group Accessions"]
            for acc in row_str.split(";")
        )

    fetch_uniprot_data(accessions)


def prefetch_all_msf_uniprot():
    """
    Fetch data for all accesions found in MSF files in MS Searched directory.
    """
    accessions = set()
    for filename in os.listdir(paths.MS_SEARCHED_DIR):
        if not os.path.splitext(filename)[1].lower() in [".msf"]:
            continue

        msf_path = os.path.join(paths.MS_SEARCHED_DIR, filename)

        with sqlite3.connect(msf_path) as conn:
            cursor = conn.cursor()

            vals = cursor.execute(
                """
                SELECT
                ProteinAnnotations.Description
                FROM
                ProteinAnnotations
                """
            )

            accessions.update(
                RE_DISCOVERER_ACCESSION.match(prot_string).group(1)
                for (prot_string,) in vals
            )

    fetch_uniprot_data(accessions)


def prefetch_accessions(psms):
    """
    Pre-fetch all UniProt information.

    Speeds up data processing by doing all queries in one go.

    Parameters
    ----------
    psms : :class:`pandas.DataFrame`
    """
    accessions = set()

    # MASCOT / Discoverer Outputs
    if "Protein Group Accessions" in psms.columns:
        accessions.update(
            acc.strip()
            for i in psms["Protein Group Accessions"]
            if isinstance(i, str)
            for acc in i.split(";")
        )

    # CAMV Outputs
    if "Accession" in psms.columns:
        accessions.update(
            acc.strip()
            for i in psms["Accession"]
            if isinstance(i, str)
            for acc in i.split("/")
        )

    fetch_uniprot_data(accessions)


def get_uniprot_data(accession):
    """
    Get UniProt data associated with a protein.

    Parameters
    ----------
    accession : str

    Returns
    -------
    dict
    """
    if accession not in UNIPROT_DATA:
        fetch_uniprot_data([accession])

    return UNIPROT_DATA[accession]
