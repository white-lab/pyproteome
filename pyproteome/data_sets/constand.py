"""
This module provides functionality for manipulating proteomics data sets.

Functionality includes merging data sets and interfacing with attributes in a
structured format.
"""

# Built-ins
from __future__ import absolute_import, division

import logging

# Core data analysis libraries
import numpy as np

import pyproteome as pyp


LOGGER = logging.getLogger("pyproteome.constand")


def constand(ds, inplace=False, n_iters=50, tol=1e-5):
    """
    Normalize channels to given levels for intra-run comparisons.

    Divides all channel values by a given level.

    Parameters
    ----------
    lvls : dict of str, float or
    ds : :class:`.DataSet`
        Mapping of channel names to normalized levels. Alternatively,
        a data set to pass to levels.get_channel_levels() or use
        pre-calculated levels from.
    inplace : bool, optional
        Modify this data set in place.

    Returns
    -------
    ds : :class:`.DataSet`
    """
    new = ds

    if not inplace:
        new = new.copy()

    LOGGER.info("Applying CONSTANd normalization.")

    channels = list(new.channels.values())

    new_channels = pyp.utils.norm(new.channels)
    # new = new.dropna()

    k = new[channels].values
    err = np.inf
    m, n = k.shape

    for ind in range(1, n_iters * 2 + 1):
        if ind % 2 == 1:
            # In the odd step the rows are fitted to match the row marginals
            # (i.e. constraints):
            #
            # K^(2t + 1) = R^(t + 1) * K ^ (2t)
            #
            # The row multipliers R^(t + 1) are computed such that the mean
            # of the reporter ion intensities equals 1/n:
            r = 1 / (n * np.nanmean(k, axis=1))

            # k^(2t + 1)
            k = np.einsum('..., ...', r, k.T).T
            err = (
                abs(np.nanmean(k, axis=0) - 1 / n)
            ).sum() / 2
        else:
            # In the even step the columns are fitted to match the column
            # marginals (i.e. constraints):
            #
            # K^(2t + 2) = K^(2t + 1) * S^(t + 1)
            #
            # The column multipliers S^(t + 1) are computed such that the
            # mean of the reporter ion intensities equals 1/n:
            s = 1 / (n * np.nanmean(k, axis=0))

            # k^(2t + 2)
            k = np.einsum('..., ...', k, s)
            err = (
                abs(
                    np.nanmean(k, axis=1) -
                    1 / np.count_nonzero(~np.isnan(k), axis=1)
                )
            ).sum() / 2

        if err < tol:
            break

    for i, (key, norm_key) in enumerate(zip(
        new.channels.values(),
        new_channels.values(),
    )):
        new.psms[norm_key] = k[:, i]
        del new.psms[key]

    new.intra_normalized = True
    new.channels = new_channels
    new.groups = new.groups.copy()

    new.update_group_changes()

    return new
