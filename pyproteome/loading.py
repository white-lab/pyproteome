"""
This module provides functionality for loading data sets.

Functionality includes loading CAMV and MASCOT / Discoverer data sets.
"""

# Built-ins
import logging

# Core data analysis libraries
import numpy as np

from . import camv, discoverer, protein, sequence, utils


LOGGER = logging.getLogger("pyproteome.loading")


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

        if not seq:
            return 0, False

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
