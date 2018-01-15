"""
Plot calculated levels of a given sequence across channels or groups.
"""

from __future__ import division

# Built-ins
from collections import OrderedDict
import itertools
import logging
import os
import re

# Core data analysis libraries
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

import pyproteome


LOGGER = logging.getLogger("pyproteome.plot")


def plot(
    data, f=None,
    title=None,
    figsize=(12, 8),
):
    """
    Plot the levels of a sequence across multiple channels.

    Parameters
    ----------
    data : :class:`DataSet<pyproteome.data_sets.DataSet>`
    sequence : str or :class:`Sequence<pyproteome.sequence.Sequence>`
    title : str, optional
    figsize : tuple of int, int
    """
    channel_names = [
        channel_name
        for group in data.groups.values()
        for channel_name in group
        if channel_name in data.channels
    ]
    channel_names = list(OrderedDict.fromkeys(channel_names))
    channels = [
        data.channels[channel_name]
        for channel_name in channel_names
    ]

    if f:
        data = data.filter(f)

    figures = []

    for _, row in data.psms.iterrows():
        seq = str(row["Sequence"])

        values = row[channels]
        mask = ~pd.isnull(row[channels])

        names = pd.Series(channel_names, index=values.index)[mask]
        values = values[mask]

        fig, ax = plt.subplots(figsize=figsize)

        indices = np.arange(len(values))
        bar_width = .35
        ax.bar(indices, values.as_matrix(), bar_width)
        ax.set_xticks(indices)
        ax.set_xticklabels(
            names,
            fontsize=20,
            rotation=45,
            horizontalalignment="right",
        )

        ax.set_title(seq + (" - {}".format(title) if title else ""))

        ax.set_ylabel(
            "Cummulative Intensity" +
            (" (Normalized)" if data.intra_normalized else "")
        )

        ax.title.set_fontsize(28)
        ax.yaxis.label.set_fontsize(20)

        fig.savefig(
            os.path.join(
                data.name,
                re.sub(
                    "[?/]",
                    "_",
                    "{}-{}.png".format(
                        " / ".join(row["Proteins"].genes),
                        data.name,
                    ),
                ),
            ),
            bbox_inches="tight",
            dpi=pyproteome.DEFAULT_DPI,
            transparent=True,
        )

        figures.append((fig, ax))

    return figures


def plot_group(
    data,
    f=None,
    cmp_groups=None,
    figsize=(12, 8),
):
    """
    Plot the levels of a sequence across each group.

    Parameters
    ----------
    data : :class:`DataSet<pyproteome.data_sets.DataSet>`
    sequences : list of str
    """
    if cmp_groups is None:
        cmp_groups = [list(data.groups.keys())]

    channels = [
        [
            [
                data.channels[channel_name]
                for channel_name in data.groups[label]
                if channel_name in data.channels
            ]
            for label in groups
        ]
        for groups in cmp_groups
    ]

    if f:
        data = data.filter(f)

    figures = []

    for _, row in data.psms.iterrows():
        names, values = [], []

        for groups in cmp_groups:
            group_vals = [
                data[[
                    data.channels[name]
                    for name in data.groups[group]
                    if name in data.channels
                ]]
                for group in groups
            ]

            groups, group_vals = zip(*[
                (group, vals.mean(axis=0))
                for group, vals in zip(groups, group_vals)
                if vals.shape[1] > 0
            ])

            # Check normalization channel is not null or all other channels
            # are not null
            if len(cmp_groups) > 1 and (
                pd.isnull(group_vals[0]).all() or all(
                    pd.isnull(vals) for vals in group_vals[1:]
                )
            ):
                continue

            groups, group_vals = zip(*[
                (group, vals[~pd.isnull(vals)])
                for group, vals in zip(groups, group_vals)
                if not pd.isnull(vals).all()
            ])

            if not groups:
                continue

            group_vals = [
                vals / group_vals[0].mean()
                for vals in group_vals
            ]

            names.append(groups)
            values.append(group_vals)

        labels = [
            name
            for group in names
            for name in group
        ]

        means = [
            i.mean()
            for vals in values
            for i in vals
        ]

        errs = [
            i.std()
            for vals in values
            for i in vals
        ]

        fig, ax = plt.subplots(figsize=figsize)

        indices = np.arange(len(means))
        width = .75
        ind = 0
        ax.bar(
            width * ind - width / 2 + indices,
            means,
            width=width,
            yerr=errs,
            # error_kw=dict(ecolor='k', lw=.25, capsize=.25, capthick=.25),
        )

        ax.set_ylabel(
            "{} Signal".format(
                "Relative" if cmp_groups else "Cumulative",
            ),
            fontsize=20,
        )
        ax.ticklabel_format(style="plain")

        for label in ax.get_yticklabels():
            label.set_fontsize(20)

        ax.set_xticks(indices - width / 2)
        ax.set_xticklabels(labels, fontsize=20)

        title = "{}".format(
            # " / ".join(sequences),
            " / ".join(
                sorted(set(row["Proteins"].genes))
            )
        )[:50]
        ax.set_title(title, fontsize=32)
        ax.xaxis.grid(False)

        # if not means[0].any():
        #     return f

        # if len(means[0]) > 1:
        #     pep_seqs = row["Sequence"]
        #     ax.legend([
        #         "{} ({} - {})".format(
        #             str(seq),
        #             seq.protein_matches[0].rel_pos,
        #             seq.protein_matches[0].rel_pos + len(seq.pep_seq),
        #         )
        #         for seq in pep_seqs
        #     ])
        #     return f, ax

        y_max = max(means) + max(errs)

        def stars(p):
            if p < 0.0001:
                return "****"
            elif (p < 0.001):
                return "***"
            elif (p < 0.01):
                return "**"
            elif (p < 0.05):
                return "*"
            else:
                return "-"

        offset = 0

        for (
            (index_a, values_a, label_a), (index_b, values_b, label_b),
        ) in itertools.combinations(zip(indices, values, labels), 2):
            if cmp_groups and not any(
                label_a in group and
                label_b in group
                for group in cmp_groups
            ):
                continue

            print(values_a, values_b)
            pval = ttest_ind(values_a, values_b).pvalue

            if pval < 0.05:
                ax.set_ylim(
                    ymax=max([
                        ax.get_ylim()[1],
                        y_max * (1 + offset * 0.1 + 0.15),
                    ]),
                )
                ax.annotate(
                    "",
                    xy=(
                        index_a + width * 1.5,
                        y_max * (1.05 + offset * 0.1)
                    ),
                    xytext=(
                        index_b + width * 1.5,
                        y_max * (1.05 + offset * 0.1)
                    ),
                    xycoords='data',
                    textcoords='data',
                    arrowprops=dict(
                        arrowstyle="-",
                        ec='#000000',
                    ),
                )
                ax.text(
                    x=np.mean([index_a, index_b]) + width * 1.5,
                    y=y_max * (1.07 + offset * 0.1),
                    s=stars(pval),
                    horizontalalignment='center',
                    verticalalignment='center',
                )

                offset += 1

            fig.savefig(
                os.path.join(
                    data.name,
                    re.sub(
                        "[?/]",
                        "_",
                        "{}-{}-between.png".format(
                            " / ".join(row["Proteins"].genes),
                            data.name,
                        ),
                    ),
                ),
                bbox_inches="tight",
                dpi=pyproteome.DEFAULT_DPI,
                transparent=True,
            )
            figures.append((fig, ax))

    return figures


def plot_all(
    data,
    f=None,
    figsize=(16, 8),
    individual=True,
    between=True,
    cmp_groups=None,
):
    try:
        os.makedirs(data.name)
    except:
        pass

    figures = []

    if individual:
        figures += plot(
            data,
            f=f,
            figsize=figsize,
        )

    if between:
        figures += plot_group(
            data,
            f=f,
            cmp_groups=cmp_groups,
        )

    return figures


def plot_together(
    data, only=True, **kwargs
):
    try:
        os.makedirs(data.name)
    except:
        pass

    prot = " / ".join(data.genes)
    figures = plot_group(
        data,
        **kwargs
    )

    return figures
