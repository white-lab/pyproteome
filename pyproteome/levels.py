'''
This module provides functionality for normalizing protein data.

Levels can be extracted from supernatant or phosphotyrosine runs using median
or mean peptide levels across multiple channels.
'''

from __future__ import absolute_import, division

# Built-ins
from collections import OrderedDict
import logging
import warnings

# Core data analysis libraries
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats

LOGGER = logging.getLogger('pyproteome.levels')
WARN_PEP_CUTOFF = 50
REL_CUTOFF = 5


def kde_max(points):
    '''
    Estimate the center of a quantification channel by fitting a gaussian
    KDE function and finding its maximum.

    Parameters
    ----------
    points : list of float

    Returns
    -------
    float
    '''
    points = points[~np.isnan(points)]
    gaus = stats.kde.gaussian_kde(points, bw_method='silverman')
    # .25 - .75 1000 slices
    # x = np.arange(0, 10, .01)
    x = np.linspace(np.quantile(points, .15), np.quantile(points, .85), 1000)
    y = np.array(gaus.pdf(x))
    return x[y == y.max()][0]


def get_channel_levels(
    data,
    norm_channels=None,
    method='median',
    cols=2,
):
    '''
    Calculate channel normalization levels. This value is calculated by
    selecting the peak of Gaussian KDE distribution fitted to channel ratio
    values.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.DataSet`
    norm_channels : list of str, optional
        Sample names of channels to use for normalization.
    method : str, optional
        Normalize to the 'mean' or 'median' of each row.
    cols : int, optional
        Number of columns used when displaying KDE distributions.

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
    channel_levels : dict of str, float
    '''
    if norm_channels is None:
        norm_channels = list(data.channels.keys())

    channels = [data.channels[i] for i in norm_channels]
    channel_levels = OrderedDict()

    rows = int(np.ceil(len(data.channels) / cols))
    scale = 4
    f, axes = plt.subplots(
        rows, cols,
        sharex=True,
        sharey=True,
        figsize=(scale * cols, scale * rows),
    )
    axes = [i for j in axes for i in j]
    ax_iter = iter(axes)

    if method in ['mean']:
        norm = data.psms[channels].mean(axis=1)
    elif method in ['median']:
        norm = data.psms[channels].median(axis=1)
    else:
        raise Exception('Unknown normalization method: {}'.format(method))

    for col_name, col in zip(norm_channels, channels):
        points = (data.psms[col] / norm).dropna()

        # Remove peptides that change more than 25x
        points = points[
            (points >= 1 / REL_CUTOFF) &
            (points <= REL_CUTOFF)
        ]

        if points.shape[0] < WARN_PEP_CUTOFF:
            LOGGER.warning(
                (
                    '{}: Too few peptides for normalization, '
                    'quantification may be inaccurate '
                    ' ({} peptides for {}: {})'
                ).format(data.name, points.shape[0], col_name, col)
            )

        if points.shape[0] < 1:
            channel_levels[col] = 1
            continue
        else:
            # Fit a guassian and find its maximum
            max_x = kde_max(points)

        channel_levels[col] = max_x

        ax = next(ax_iter)

        # seaborn==0.9.0 throws a scipy.stats warning
        with warnings.catch_warnings():
            warnings.filterwarnings(
                'ignore',
                '',
                FutureWarning,
            )
            sns.distplot(
                points.apply(np.log2),
                bins=25,
                ax=ax,
            )

        ax.set_xlim(left=-2, right=2)

        ax.set_title(
            '{} ({})'.format(col_name, col)
            if isinstance(data.channels, dict) else
            col,
        )

        txt = 'center = {:.2f}\n$\\sigma$ = {:.2f}'.format(
            max_x,
            points.std(ddof=1),
        )
        lines = [
            ax.axvline(np.log2(points.median()), color='g', zorder=2, linestyle=':'),
            ax.axvline(np.log2(points.mean()), color='r', zorder=3, linestyle='-.'),
            ax.axvline(np.log2(max_x), color='c', zorder=4, linestyle='--'),
            ax.axvline(np.log2(1), color='k', linestyle='-', zorder=0, lw=2),
        ]
        if col_name is norm_channels[0]:
            ax.legend(lines[:-1], ['kde max', 'mean', 'median'])

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
                facecolor='white',
                zorder=1,
                edgecolor='black',
                boxstyle='round',
            )
        )

    for ax in ax_iter:
        ax.axis('off')

    f.suptitle(
        '{}'.format(data.name),
        fontsize=16,
    )

    return f, channel_levels
