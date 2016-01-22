"""
This module provides functionality for processing data sets.

Functionality includes normalizing data, filtering bad peptide hits, and
calculating changes between groups.
"""
# Built-ins
import logging

# Core data analysis libraries
import pandas as pd

from . import utils, math


LOGGER = logging.getLogger("pyproteome.processing")


def filter_mascot_matches(
    psms, channel_names,
    confidence="High", ion_score_cutoff=20,
):
    """
    Filter a list of MASCOT matches.

    This filter only keeps those above a given quality threshold.

    Parameters
    ----------
    psms : pandas.DataFrame
    channel_names : list of str
    confidence : str, optional
    ion_score_cutoff : int or None, optional

    Returns
    -------
    pandas.DataFrame
    """
    if confidence and "Confidence Level" in psms.columns:
        psms = psms[psms["Confidence Level"] == "High"]

    if ion_score_cutoff is not None:
        psms = psms[psms["IonScore"] >= ion_score_cutoff]

    psms = psms.dropna(
        axis=0,
        how="any",
        subset=channel_names,
    )

    return psms


def preprocess_mascot_matches(psms, channel_names):
    """
    Process a table of MASCOT matches.

    Turns proteins, sequences, and modifications into internal python objects.
    Also groups together duplicate hits, summing their channel levels together.

    Parameters
    ----------
    psms : pandas.DataFrame
    channel_names : list of str

    Returns
    -------
    pandas.DataFrame
    """
    grouped = psms.groupby(
        [
            "Proteins",
            "Sequence",
            "Modifications",
        ],
        sort=False,
        as_index=False,
    )
    agg_dict = dict(
        (i, sum)
        for i in channel_names
    )
    processed = grouped.agg(agg_dict)
    return processed


def _normalize_matches(psms, levels, group_a=None, group_b=None):
    """
    Normalize a list of peptides to given levels.

    Levels should be calculated using a function defined within
    pyproteome.levels.

    Parameters
    ----------
    psms : pandas.DataFrame
    levels : dict of str, float
    group_a : list of str, optional
    group_b : list of str, optional

    Returns
    -------
    pandas.DataFrame

    See Also
    --------
    pyproteome.levels
    """
    for col, level in levels.items():
        psms[utils.norm(col)] = psms[col] / level

    if group_a and group_b:
        psms = _update_snr_change(psms, group_a, group_b)

    return psms


def normalize_psms(
    psms, channel_levels,
    filter_psms=True,
    group_a=None, group_b=None,
):
    """
    Normalize a list of peptides to given channel levels.

    Channel levels should be calculated using a function defined within
    pyproteome.levels.

    Parameters
    ----------
    psms : pandas.DataFrame
    levels : dict of str, float
    group_a : list of str, optional
    group_b : list of str, optional

    Returns
    -------
    pandas.DataFrame

    See Also
    --------
    pyproteome.levels
    """
    channel_names = list(channel_levels.keys())

    # Filter Unreliable PSMs
    if filter_psms:
        filtered_psms = filter_mascot_matches(
            psms,
            channel_names,
            confidence="High",
            ion_score_cutoff=None,
        )
    else:
        filtered_psms = filter_mascot_matches(
            psms,
            channel_names,
            confidence=None,
            ion_score_cutoff=None,
        )

    # Pre-process Matches to remove duplicate rows
    processed_psms = preprocess_mascot_matches(filtered_psms, channel_names)

    # Normalize Channel Levels
    norm_psms = _normalize_matches(
        processed_psms, channel_levels,
        group_a=group_a, group_b=group_b
    )

    print(norm_psms.shape)
    return norm_psms


def _update_snr_change(psms, group_a, group_b, normalize=True):
    """
    Update a table's SNR, Fold-Change, and p-value columns.

    Values are calculated based on changes between group_a and group_b.

    Parameters
    ----------
    psms : pandas.DataFrame
    group_a : list of str
    group_b : list of str
    normalize : bool, optional

    Returns
    -------
    pandas.DataFrame
    """
    if normalize:
        group_a = utils.norm(group_a)
        group_b = utils.norm(group_b)

    psms = psms.copy()

    psms["SNR"] = pd.Series(
        [
            math.snr(row[group_a], row[group_b])
            for index, row in psms.iterrows()
        ]
    )
    psms["Fold Change"] = pd.Series(
        [
            math.fold_change(row[group_a], row[group_b])
            for index, row in psms.iterrows()
        ]
    )
    psms["p-value"] = pd.Series(
        [
            math.p_value(row[group_a], row[group_b])
            for index, row in psms.iterrows()
        ]
    )
    return psms


def merge_tables(psms_lst, group_a, group_b):
    """
    Merge together two datasets.

    Sums peptide levels in cases of duplicates.

    Parameters
    ----------
    psms_lst : list of pandas.DataFrame
    group_a : list of str
    group_b : list of str

    Returns
    -------
    pandas.DataFrame
    """
    a_norm = utils.norm(group_a)
    b_norm = utils.norm(group_b)

    agg_dict = dict((i, sum) for i in group_a + group_b + a_norm + b_norm)

    output = pd.concat(psms_lst).groupby(
        [
            "Proteins",
            "Sequence",
            "Modifications",
        ],
        sort=False,
        as_index=False,
    ).agg(agg_dict)

    output = _update_snr_change(output, group_a, group_b)

    return output
