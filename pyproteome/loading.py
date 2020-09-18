'''
This module provides functionality for loading data sets.

Functionality includes loading CAMV and Proteome Discoverer data sets.
'''

# Built-ins
import logging

# Core data analysis libraries
import numpy as np

from . import camv, discoverer


LOGGER = logging.getLogger('pyproteome.loading')


def _calculate_rejected(psms, accepted, maybed, rejected):
    if rejected is None:
        return psms

    LOGGER.info('Filtering out rejected scans.')

    # Remove any peptides that match the scan number and sequence
    # in the rejected list.
    reject_mask = np.zeros(psms.shape[0], dtype=bool)
    validations = psms['Validated'].copy()

    for index, row in psms.iterrows():
        # Check if this specific sequence and scan was rejected
        hit = np.logical_and(
            # Assuming Scan always == Last Scan
            rejected['Scan'] == row['Scan'],
            rejected['Sequence'] == row['Sequence'],
        )

        if hit.any():
            reject_mask[index] = True
            continue

        if accepted is not None:
            if np.logical_and(
                accepted['Scan'] == row['Scan'],
                accepted['Sequence'] == row['Sequence'],
            ).any():
                validations[index] = True
                continue

        if maybed is not None:
            if np.logical_and(
                maybed['Scan'] == row['Scan'],
                maybed['Sequence'] == row['Sequence'],
            ).any():
                continue

        # Check if this scan was rejected and no sequences were accepted
        hit = (rejected['Scan'] == row['Scan']).any()
        if not hit:
            continue

        reject_mask[index] = True

    psms['Validated'] = validations
    psms = psms[~reject_mask].reset_index(drop=True)

    return psms


def _calculate_accepted(psms, accepted):
    if accepted is None:
        return psms

    LOGGER.info('Filtering out non-accepted scans.')

    reject_mask = np.zeros(psms.shape[0], dtype=bool)
    validations = psms['Validated'].copy()

    for index, row in psms.iterrows():
        # Reject hits where the scan number is the same but the sequence
        # is different.
        hit = np.logical_and(
            accepted['Scan'] == row['Scan'],
            accepted['Sequence'] != row['Sequence'],
        )
        if hit.any():
            reject_mask[index] = True

        hit = np.logical_and(
            accepted['Scan'] == row['Scan'],
            accepted['Sequence'] != row['Sequence'],
        )
        if hit.any():
            validations[index] = True

    psms['Validated'] = validations
    psms = psms[~reject_mask].reset_index(drop=True)

    return psms


def load_psms(basename, pick_best_psm=True):
    '''
    Load a list of peptide-spectrum matches (PSMs) from a .msf file produced by
    Proteome Discoverer.

    Parameters
    ----------
    basename : str
        Base name of the data set (i.e. 'CK-H1-pY' for 'CK-H1-pY.msf').
    pick_best_psm : bool, optional
        Select the best scoring PSM for a given scan, otherwise load all PSMs.

    Returns
    -------
    psms : :class:`pandas.DataFrame`
    '''
    # The load CAMV data to clear unwanted hits if available.
    accepted, maybed, rejected = camv.load_camv_validation(basename)
    lst = (accepted, maybed, rejected)

    psms, species = discoverer.read_discoverer_msf(
        basename,
        pick_best_psm=(
            pick_best_psm and
            all([not i for i in lst])
        ),
    )

    psms = _calculate_rejected(psms, accepted, maybed, rejected)
    psms = _calculate_accepted(psms, accepted)

    return psms, species, lst
