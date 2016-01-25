"""This module provides functions for mathematical calculations."""

# Core data analysis libraries
import numpy as np


def log_cum_fold_change(vals):
    """
    Calculate the cumulative fold change (in base-2) of values.

    Fold-change is normalized to the first element in the array.

    Parameters
    ----------
    vals : numpy.array

    Returns
    -------
    float
    """
    return sum(
        abs(np.log2(i / (vals[0])))
        for i in vals[1:]
        if i != 0
    )


def snr(data1, data2):
    """
    Calculate the signal-to-noise ratio between two groups of data.

    Parameters
    ----------
    data1 : numpy.array of float
    data2 : numpy.array of float

    Returns
    -------
    float
    """
    return (data1.mean() - data2.mean()) \
        / (data1.std(ddof=1) + data2.std(ddof=1))
