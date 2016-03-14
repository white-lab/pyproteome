"""
This module provides functionality for fetching protein data from UniProt.

Caches fetched protein data for faster re-use.
"""

# Built-ins
import os
import re

# Core data analysis libraries
import pandas as pd

import uniprot

from . import paths


RE_ACCESSION = re.compile("\[([A-Za-z0-9]+_[A-Z]+)\]")

UNIPROT_DATA = {}
MISS_COUNT = {}


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

    UNIPROT_DATA.update(
        uniprot.get_metadata_with_some_seqid_conversions(
            accessions,
            cache_dir=os.path.abspath('.uniprot_cache'),
        )
    )


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
