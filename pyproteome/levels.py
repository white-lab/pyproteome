"""
This module provides functionality for normalizing protein data.

Levels can be extracted from supernatant or phosphotyrosine runs using median
or mean peptide levels across multiple channels.
"""

from __future__ import division

# Built-ins
import os
from collections import OrderedDict

# Core data analysis libraries
from matplotlib import pyplot as plt
import numpy as np

from . import utils


def get_channel_levels(
    data,
    folder_name=None, file_name=None,
):
    """
    Calculate channel levels using median channel levels.

    Should be used on sample supernatants.

    Parameters
    ----------
    data : :class:`DataSet<pyproteome.data_sets.DataSet>`
    folder_name : str, optional
    file_name : str, optional
    """
    if folder_name is None:
        folder_name = data.name

    if not file_name:
        file_name = "channel_levels.png"

    utils.make_folder(folder_name)

    if folder_name:
        file_name = os.path.join(folder_name, file_name)

    channel_names = list(data.channels.keys())
    channels = list(data.channels.values())
    channel_levels = OrderedDict()
    base = channels[0]
    channel_levels[base] = 1

    f, axes = plt.subplots(
        int(np.ceil(len(data.channels) / 2)), 2,
        sharex=True,
        figsize=(12, 12),
    )
    axes = [i for j in axes for i in j]

    for ax, col_name, col in zip(list(axes), channel_names[1:], channels[1:]):
        points = (data.psms[col] / data.psms[base]).dropna().as_matrix()

        # Filter ratios > 30, those are likely an error in quantification and
        # may throw off histogram binning
        points = points[points < 30]
        med = np.median(points)
        channel_levels[col] = med

        ax.hist(points, bins=40)
        ax.set_title(
            "{}: median: {:.2f}, $\sigma$ = {:.2f}".format(
                "{} ({})".format(col_name, col)
                if isinstance(data.channels, dict) else
                col,
                med,
                points.std(ddof=1),
            )
        )
        ax.axvline(med, color='k', linestyle='--')

    for ax in axes[len(data.channels) - 1:]:
        ax.set_axis_off()

    f.show()

    if file_name:
        f.savefig(file_name)

    return channel_levels


def get_average_phospho_levels(data):
    """
    Calculate channel levels using mean channel levels.

    Parameters
    ----------
    data : :class:`DataSet<pyproteome.data_sets.DataSet>`

    Returns
    -------
    OrderedDict
    """
    base = list(data.channels.keys())[0]
    levels = data.psms[list(data.channels)].mean()

    return OrderedDict(
        (name, levels[name] / levels[base])
        for name in data.channels
    )
