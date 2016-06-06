"""
This module provides functionality for loading data sets.

Functionality includes loading CAMV and MASCOT / Discoverer data sets.
"""
# Built-ins
import logging
import os
import re

# Core data analysis libraries
import numpy as np
import pandas as pd

from . import (
    camv, discoverer, fetch_data, protein, sequence, modification, paths, utils
)


LOGGER = logging.getLogger("pyproteome.loading")
RE_PROTEIN = re.compile("([A-Za-z0-9\(\)\[\]\\/\',\. \-\+]+) OS=")
RE_MODIFICATION = re.compile("(N-Term|C-Term|([A-Z])([0-9]*))\((.*)\)")


def _filter_unassigned_rows(psms):
    """
    Remove rows from psms with unassigned peptides / proteins.

    Parameters
    ----------
    psms : :class:`pandas.DataFrame`

    Returns
    -------
    :class:`pandas.DataFrame`
    """
    return psms.dropna(
        axis=0,
        how="any",
        subset=[
            "Protein Descriptions",
            "Protein Group Accessions",
            "Sequence"
        ],
    )


def _extract_protein(prot_string):
    match = fetch_data.RE_ACCESSION.search(prot_string)
    if not match:
        raise Exception(
            "Unable to find accession in \"{}\"".format(prot_string)
        )
    return protein.Protein(accession=match.group(0))


def _extract_proteins_from_accessions(prots_string):
    """
    Extract a list of proteins from a string containing protein accessions.

    Parameters
    ----------
    prots_string : str

    Returns
    -------
    :class:`Proteins<pyproteome.protein.Proteins`
    """
    return protein.Proteins(
        proteins=[
            protein.Protein(
                accession=accession.strip(),
            )
            for accession in prots_string.split(";")
        ],
    )


def _extract_proteins_from_description(prots_string):
    """
    Extract a list of proteins from a string containing protein descriptions.

    Parameters
    ----------
    prots_string : str

    Returns
    -------
    :class:`Proteins<pyproteome.protein.Proteins`
    """
    return protein.Proteins(
        proteins=[
            protein.Protein(
                accession=accession,
            )
            for accession in fetch_data.RE_ACCESSION.findall(prots_string)
        ],
    )


def extract_sequence(proteins, sequence_string):
    """
    Extract a Sequence object from a list of proteins and sequence string.

    Does not set the Sequence.modifications attribute.

    Parameters
    ----------
    proteins : list of :class:`Protein<pyproteome.protein.Protein>`
    sequence_string : str

    Returns
    -------
    list of :class:`Sequence<pyproteome.sequence.Sequence>`
    """
    prot_matches = []

    # Skip peptides with no protein matches
    if not isinstance(proteins, protein.Proteins):
        proteins = []

    def _get_rel_pos(protein, pep_seq):
        seq = protein.full_sequence
        pep_pos = seq.find(pep_seq)
        exact = True

        if pep_pos < 0:
            pep_pos = utils.fuzzy_find(pep_seq, seq)
            exact = False

        return pep_pos, exact

    for prot in proteins:
        rel_pos, exact = _get_rel_pos(prot, sequence_string.upper())

        prot_matches.append(
            sequence.ProteinMatch(
                protein=prot,
                rel_pos=rel_pos,
                exact=exact,
            )
        )

    return sequence.Sequence(
        pep_seq=sequence_string,
        protein_matches=prot_matches,
    )


def _extract_modification(seq, mod_string):
    """
    Extract a single modification from mod_string.

    Parameters
    ----------
    seq : :class:`Sequence<pyproteome.sequence.Sequence>`
    mod_string : str

    Returns
    -------
    :class:`Modification<pyproteome.modification.Modification>`
    """
    def _get_pos(match):
        if match.group(1) == "N-Term":
            return 0, True, False
        elif match.group(1) == "C-Term":
            return len(seq.pep_seq) - 1, False, True
        else:
            return int(match.group(3)) - 1, False, False

    match = RE_MODIFICATION.match(mod_string)
    pos, nterm, cterm = _get_pos(match)
    letter, mod = match.group(2), match.group(4)

    if not nterm and not cterm:
        assert seq.pep_seq[pos].upper() == letter

    return modification.Modification(
        rel_pos=pos,
        mod_type=mod,
        nterm=nterm,
        cterm=cterm,
        sequence=seq,
    )


def _extract_modifications(sequence, mods_string):
    """
    Extract a structured list of modifications from mod_string.

    Parameters
    ----------
    sequence : :class:`Sequence<pyproteome.sequence.Sequence>`
    mod_string : str

    Returns
    -------
    :class:`Modifications<pyproteome.modification.Modifications>`
    """
    if isinstance(mods_string, float):
        mods_string = ""

    return modification.Modifications(
        mods=[
            _extract_modification(sequence, i.strip())
            for i in mods_string.split(";")
            if i
        ],
    )


def read_table_delimited(basename):
    psms_path = os.path.join(
        paths.MS_SEARCHED_DIR,
        basename + "_psms.txt",
    )

    LOGGER.info(
        "Loading MASCOT peptides from \"{}\"".format(
            os.path.basename(psms_path),
        )
    )

    psms = pd.read_table(psms_path)
    psms = _filter_unassigned_rows(psms)

    # Pre-fetch UniProt data to speed up later queries
    fetch_data.prefetch_accessions(psms)

    psms["Proteins"] = pd.Series(
        [
            _extract_proteins_from_accessions(
                row["Protein Group Accessions"]
            )
            for index, row in psms.iterrows()
        ],
        index=psms.index,
    )
    psms["Sequence"] = pd.Series(
        [
            extract_sequence(row["Proteins"], row["Sequence"])
            for index, row in psms.iterrows()
        ],
        index=psms.index,
    )
    psms["Modifications"] = pd.Series(
        [
            _extract_modifications(row["Sequence"], row["Modifications"])
            for index, row in psms.iterrows()
        ],
        index=psms.index,
    )

    psms.reset_index(inplace=True, drop=True)

    # Finally close the reference loop between sequences and modifications
    for index, row in psms.iterrows():
        row["Sequence"].modifications = row["Modifications"]

    return psms


def _calculate_rejected(psms, accepted, maybed, rejected):
    if rejected is None:
        return psms

    LOGGER.info("Filtering out rejected scans.")

    # Remove any peptides that match the scan number and sequence
    # in the rejected list.
    reject_mask = np.zeros(psms.shape[0], dtype=bool)
    validations = psms["Validated"].copy()

    for index, row in psms.iterrows():
        # Check if this specific sequence and scan was rejected
        hit = np.logical_and(
            # Assuming First Scan always == Last Scan
            rejected["Scan"] == row["First Scan"],
            rejected["Sequence"] == row["Sequence"],
        )

        if hit.any():
            reject_mask[index] = True
            continue

        if accepted is not None:
            if np.logical_and(
                accepted["Scan"] == row["First Scan"],
                accepted["Sequence"] == row["Sequence"],
            ).any():
                validations[index] = True
                continue

        if maybed is not None:
            if np.logical_and(
                maybed["Scan"] == row["First Scan"],
                maybed["Sequence"] == row["Sequence"],
            ).any():
                continue

        # Check if this scan was rejected and no sequences were accepted
        hit = (rejected["Scan"] == row["First Scan"]).any()
        if not hit:
            continue

        reject_mask[index] = True

    psms["Validated"] = validations
    psms = psms[~reject_mask].reset_index(drop=True)

    return psms


def _calculate_accepted(psms, accepted):
    if accepted is None:
        return psms

    LOGGER.info("Filtering out non-accepted scans.")

    reject_mask = np.zeros(psms.shape[0], dtype=bool)
    validations = psms["Validated"].copy()

    for index, row in psms.iterrows():
        # Reject hits where the scan number is the same but the sequence
        # is different.
        hit = np.logical_and(
            accepted["Scan"] == row["First Scan"],
            accepted["Sequence"] != row["Sequence"],
        )
        if hit.any():
            reject_mask[index] = True

        hit = np.logical_and(
            accepted["Scan"] == row["First Scan"],
            accepted["Sequence"] != row["Sequence"],
        )
        if hit.any():
            validations[index] = True

    psms["Validated"] = validations
    psms = psms[~reject_mask].reset_index(drop=True)

    return psms


def load_mascot_psms(basename, camv_slices=None, msf=True):
    """
    Load a list of sequences from a MSF file produced by MASCOT / Discoverer.

    Parameters
    ----------
    basenme : str
    camv_slices : int, optional

    Returns
    -------
    psms : :class:`pandas.DataFrame`
    scan_lists : dict of str, list of int
    msf : bool, optional
        Use Discoverer .msf files. Otherwise look for tab-delimited files.
    """
    if msf:
        psms = discoverer.read_discoverer_msf(basename)
        psms = _filter_unassigned_rows(psms)
    else:
        psms = read_table_delimited(basename)

    # Output the phosphorylation scan list for CAMV
    psms, scan_lists = camv.output_scan_list(
        psms,
        basename=basename,
        letter_mod_types=[(None, "Phospho")],
        scan_sets=camv_slices,
    )

    psms["Validated"] = False

    # The load CAMV data to clear unwanted hits if available.
    accepted, maybed, rejected = camv.load_camv_validation(basename)

    psms = _calculate_rejected(psms, accepted, maybed, rejected)
    psms = _calculate_accepted(psms, accepted)

    return psms, scan_lists


def load_validated_psms(filename):
    """
    Wrap load_camv_validation, returning only accepted peptides.

    Parameters
    ----------
    filename : str

    Returns
    -------
    :class:`pandas.DataFrame`
    """
    accepted, maybed, rejected = camv.load_camv_validation(filename)

    accepted["Validated"] = True
    return accepted
