"""
This module provides functionality for normalizing protein data.

Levels can be extracted from supernatant or phosphotyrosine runs using median
or mean peptide levels across multiple channels.
"""

from __future__ import absolute_import, division

# Built-ins
from collections import OrderedDict
import logging
import os
import warnings

# Core data analysis libraries
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats

import pyproteome as pyp

LOGGER = logging.getLogger("pyproteome.levels")
WARN_PEP_CUTOFF = 50


def get_channel_levels(
    data,
    folder_name=None,
    file_name=None,
    cols=2,
):
    """
    Calculate channel normalization levels. This value is calculated by
    selecting the peak of Gaussian KDE distribution fitted to channel ratio
    values.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.DataSet`
    folder_name : str, optional
    file_name : str, optional
    cols : int, optional
        Number of columns used when displaying KDE distributions.

    Returns
    -------
    dict of str, float
    """
    if not file_name:
        file_name = "channel_levels.png"

    folder_name = pyp.utils.make_folder(
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
        sharey=True,
        figsize=(3 * cols, 3 * rows),
    )
    axes = [i for j in axes for i in j]
    ax_iter = iter(axes)

    means = data.psms[channels].mean(axis=1)
    # return {key: 1 for key in data.channels.values()}

    for col_name, col in zip(channel_names, channels):
        points = (data.psms[col] / means).dropna()

        if points.shape[0] < WARN_PEP_CUTOFF:
            LOGGER.warning(
                (
                    "{}: Too few peptides for normalization, "
                    "quantification may be inaccurate "
                    " ({} peptides for {}: {})"
                ).format(data.name, points.shape[0], col_name, col)
            )

        if points.shape[0] < 1:
            channel_levels[col] = 1
            continue
        else:
            # Fit a guassian and find its maximum
            gaus = stats.kde.gaussian_kde(points)
            x = np.arange(0, 10, .01)
            y = np.array(gaus.pdf(x))
            med = x[y == y.max()][0]

        channel_levels[col] = med

        ax = next(ax_iter)

        # seaborn==0.9.0 throws a scipy.stats warning
        with warnings.catch_warnings():
            warnings.filterwarnings(
                'ignore',
                '',
                FutureWarning,
            )
            sns.distplot(
                points,
                bins=25,
                ax=ax,
            )

        ax.set_title(
            "{} ({})".format(col_name, col)
            if isinstance(data.channels, dict) else
            col,
        )

        txt = "center = {:.2f}\n$\\sigma$ = {:.2f}".format(
            med,
            points.std(ddof=1),
        )
        ax.axvline(med, color='k', linestyle='--')

        ax.text(
            s=txt,
            x=ax.get_xlim()[1] * .9,
            y=1,
            color='k',
            horizontalalignment='right',
            verticalalignment='center',
        ).set_bbox(
            dict(
                # facecolor=_get_color(txt, x, y),
                alpha=1,
                linewidth=0.5,
                facecolor="white",
                zorder=1,
                edgecolor="black",
                boxstyle="round",
            )
        )

    for ax in ax_iter:
        ax.axis("off")

    f.suptitle(
        "{}".format(data.name),
        fontsize=16,
    )

    if file_name:
        f.savefig(
            os.path.join(folder_name, file_name),
            bbox_inches="tight",
            dpi=pyp.utils.DEFAULT_DPI,
            transparent=True,
        )

    return channel_levels
