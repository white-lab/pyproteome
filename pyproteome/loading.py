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


def _extract_protein(prot_string):
    match = fetch_data.RE_ACCESSION.search(prot_string)
    if not match:
        raise Exception(
            "Unable to find accession in \"{}\"".format(prot_string)
        )
    return protein.Protein(accession=match.group(0))


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


def load_mascot_psms(basename, camv_slices=None, pick_best_ptm=False):
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
    pick_best_ptm : bool, optional
    """
    # The load CAMV data to clear unwanted hits if available.
    accepted, maybed, rejected = camv.load_camv_validation(basename)
    lst = (accepted, maybed, rejected)

    psms = discoverer.read_discoverer_msf(
        basename,
        pick_best_ptm=(
            pick_best_ptm and
            all(not i for i in lst)
        ),
    )

    # Output the phosphorylation scan list for CAMV
    psms, scan_lists = camv.output_scan_list(
        psms,
        basename=basename,
        letter_mod_types=[(None, "Phospho")],
        scan_sets=camv_slices,
    )

    psms["Validated"] = False

    psms = _calculate_rejected(psms, accepted, maybed, rejected)
    psms = _calculate_accepted(psms, accepted)

    return psms, scan_lists, lst


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
