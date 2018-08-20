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
import seaborn as sns
from scipy import stats

from . import utils

LOGGER = logging.getLogger("pyproteome.levels")
WARN_PEP_CUTOFF = 100


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
    ax_iter = iter(axes)

    means = data.psms[channels].mean(axis=1)

    for col_name, col in zip(channel_names, channels):
        points = (data.psms[col] / means).dropna()

        # Filter ratios > 30, those are likely an error in quantification and
        # may throw off histogram binning
        points = points[points < 30]

        if points.shape[0] < WARN_PEP_CUTOFF:
            LOGGER.warning(
                (
                    "{}: Too few peptides for normalization, "
                    "quantification may be inaccurate "
                    " ({} peptides for {}: {})"
                ).format(data.name, points.shape[0], col_name, col)
            )

        # Fit a guassian and find its maximum
        gaus = stats.kde.gaussian_kde(points)
        x = np.arange(0, 10, .01)
        y = np.array(gaus.pdf(x))
        med = x[y == y.max()][0]

        channel_levels[col] = med

        ax = next(ax_iter)

        sns.distplot(
            points,
            ax=ax,
        )
        ax.set_title(
            r"{}: median: {:.2f}, $\sigma$ = {:.2f}".format(
                "{} ({})".format(col_name, col)
                if isinstance(data.channels, dict) else
                col,
                med,
                points.std(ddof=1),
            ),
        )
        ax.axvline(med, color='k', linestyle='--')

    for ax in ax_iter:
        ax.axis("off")

    f.suptitle(
        "{}".format(data.name),
        fontsize=20,
    )

    if file_name:
        f.savefig(
            os.path.join(folder_name, file_name),
            bbox_inches="tight",
            dpi=utils.DEFAULT_DPI,
            transparent=True,
        )

    return channel_levels
