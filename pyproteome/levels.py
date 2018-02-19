"""
This module provides functionality for normalizing protein data.

Levels can be extracted from supernatant or phosphotyrosine runs using median
or mean peptide levels across multiple channels.
"""

from __future__ import absolute_import, division

# Built-ins
import logging
import os
from collections import OrderedDict

# Core data analysis libraries
from matplotlib import pyplot as plt
import numpy as np

from . import utils

LOGGER = logging.getLogger("pyproteome.levels")
WARN_PEP_CUTOFF = 150


def get_channel_levels(
    data,
    folder_name=None,
    file_name=None,
    cols=2,
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
    if not file_name:
        file_name = "channel_levels.png"

    folder_name = utils.make_folder(
        data=data,
        folder_name=folder_name,
        sub="Normalization",
    )

    channel_names = list(data.channels.keys())
    channels = list(data.channels.values())
    channel_levels = OrderedDict()
    base = channels[0]
    channel_levels[base] = 1

    rows = int(np.ceil(len(data.channels) / cols))
    f, axes = plt.subplots(
        rows, cols,
        sharex=True,
        figsize=(2 * rows, 6 * cols),
    )
    axes = [i for j in axes for i in j]

    for ax, col_name, col in zip(list(axes), channel_names[1:], channels[1:]):
        points = (data.psms[col] / data.psms[base]).dropna().as_matrix()

        # Filter ratios > 30, those are likely an error in quantification and
        # may throw off histogram binning
        points = points[points < 30]

        if points.shape[0] < WARN_PEP_CUTOFF:
            LOGGER.warning(
                (
                    "{}: Too few peptides for normalization, "
                    "quantification may be inaccurate ({} peptides for {}: {})"
                ).format(data.name, points.shape[0], col_name, col)
            )

        med = np.median(points)
        channel_levels[col] = med

        ax.hist(points, bins=40)
        ax.set_title(
            r"{}: median: {:.2f}, $\sigma$ = {:.2f}".format(
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

    f.suptitle("{}".format(data.name))

    if file_name:
        f.savefig(
            os.path.join(folder_name, file_name),
            bbox_inches="tight",
            dpi=utils.DEFAULT_DPI,
            transparent=True,
        )

    return channel_levels
