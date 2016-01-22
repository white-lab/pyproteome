"""
This module provides functionality for normalizing protein data.

Levels can be extracted from supernatant or phosphotyrosine runs using median
or mean peptide levels across multiple channels.
"""

# Built-ins
import os
from collections import OrderedDict

# Core data analysis libraries
from matplotlib import pyplot as plt
import numpy as np

from . import utils, processing


def get_channel_levels(
    psms, channel_names,
    folder_name=None, file_name="channel_levels.png",
):
    """
    Calculate channel levels using median channel levels.

    Should be used on sample supernatants.

    Parameters
    ----------
    psms : pandas.DataFrame
    channel_names : list of str or dict of str, str
    folder_name : str, optional
    file_name : str, optional
    """
    utils.make_folder(folder_name)

    if folder_name:
        file_name = os.path.join(folder_name, file_name)

    if isinstance(channel_names, dict):
        keys = list(channel_names.keys())
    else:
        keys = channel_names

    psms = processing.filter_mascot_matches(
        psms,
        keys,
    )

    channel_levels = OrderedDict()
    channel_levels[keys[0]] = 1

    f, axes = plt.subplots(
        len(keys) // 2, 2,
        sharex=True,
        figsize=(12, 12),
    )
    axes = [i for j in axes for i in j]

    for ax, col in zip(list(axes), keys[1:]):
        data = (psms[col] / psms[keys[0]]).dropna().as_matrix()
        med = np.median(data)
        channel_levels[col] = med

        ax.hist(data, bins=40)
        ax.set_title(
            "{}: median: {:.2f}, $\sigma$ = {:.2f}".format(
                "{} ({})".format(channel_names[col], col)
                if isinstance(channel_names, dict) else
                col,
                med,
                data.std(ddof=1),
            )
        )
        ax.axvline(med, color='k', linestyle='--')

    for ax in axes[len(channel_names) - 1:]:
        ax.set_axis_off()

    f.show()

    if file_name:
        f.savefig(file_name)

    return channel_levels


def get_average_phospho_levels(psms, channel_names):
    """
    Calculate channel levels using mean channel levels.

    Parameters
    ----------
    psms : pandas.DataFrame
    channel_names : list of str or dict of str, str

    Returns
    -------
    OrderedDict
    """
    if isinstance(channel_names, dict):
        keys = list(channel_names.keys())
    else:
        keys = channel_names

    levels = psms[keys].mean()

    return OrderedDict(
        (name, levels[name] / levels[keys[0]])
        for name in keys
    )
