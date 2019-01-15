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
import seaborn as sns
from scipy.stats import ttest_ind

import pyproteome as pyp


LOGGER = logging.getLogger("pyproteome.plot")


def plot(
    data,
    title=None,
    folder_name=None,
    ax=None,
):
    """
    Plot the levels of a sequence across multiple channels.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.DataSet`
    title : str, optional
    folder_name : str, optional
    figsize : tuple of (int, int), optional

    Returns
    -------
    figs : list of :class:`matplotlib.figure.Figure`
    """
    folder_name = pyp.utils.make_folder(
        data=data,
        folder_name=folder_name,
        sub="Peptides",
    )

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

    figures = []

    for _, row in data.psms.iterrows():
        seq = str(row["Sequence"])

        values = row[channels]
        mask = ~pd.isnull(row[channels])

        names = pd.Series(channel_names, index=values.index)[mask]
        values = values[mask]
        cmp_groups = data.cmp_groups or [[
            grp
            for grp, channels in data.groups.items()
            if any([i in data.channels for i in channels])
        ]]
        values = values / (
            values[[
                data.channels[i]
                for grp in cmp_groups
                for i in data.groups[grp[0]]
                if i in values.index
            ]].median()
            if len(cmp_groups) > 0 else
            values[0]
        )

        if ax is None:
            fig, ax_i = plt.subplots(
                figsize=(len(values) / 2, 6 / 2),
            )
        else:
            ax_i = ax

        df = pd.DataFrame(
            [
                (
                    name,
                    val,
                    [
                        group_name
                        for group_name, group in data.groups.items()
                        if name in group
                    ][0],
                )
                for name, val in zip(names, values)
            ],
            columns=("name", "val", "group"),
        )
        df["val"] = df["val"].apply(np.log2)

        sns.barplot(
            x="name",
            y="val",
            hue="group",
            data=df,
            ax=ax_i,
            dodge=False,
        )
        ax_i.set_xticklabels(
            ax_i.get_xticklabels(),
            rotation=45,
            horizontalalignment="right",
        )
        ax_i.set_xlabel("")
        ax_i.get_legend().set_title("")

        ax_i.axhline(np.log2(1), linestyle=":", color="k", alpha=.5)

        mod_str = row["Modifications"].__str__(prot_index=0)

        ax_i.set_title(
            title
            if title else
            "{} ({}{})".format(
                seq,
                " / ".join(row["Proteins"].genes)[:20],
                (" " + mod_str) if mod_str else "",
            ),
        )

        ylabel = "Intensity"
        # if data.cmp_groups:
        # else:
        #     ylabel = (
        #         "Cummulative Intensity" +
        #         (" (Normalized)" if data.intra_normalized else "")
        #     )

        ax_i.set_ylabel(
            ylabel,
        )

        ax_i.get_figure().savefig(
            os.path.join(
                folder_name,
                re.sub(
                    "[?/]",
                    "_",
                    "{} - {} - all.png".format(
                        "+".join(row["Proteins"].genes)[:100],
                        row["Sequence"],
                    ),
                ),
            ),
            bbox_inches="tight",
            dpi=pyp.DEFAULT_DPI,
            transparent=True,
        )

        figures.append((ax_i.get_figure(), ax_i))

    return figures


def plot_group(
    data,
    cmp_groups=None,
    title=None,
    folder_name=None,
    ax=None,
):
    """
    Plot the levels of a sequence across each group.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.DataSet`
    cmp_groups : list of tuple, optional
    title : str, optional
    folder_name : str, optional
    figsize : tuple of (int, int), optional

    Returns
    -------
    figs : list of :class:`matplotlib.figure.Figure`
    """
    folder_name = pyp.utils.make_folder(
        data=data,
        folder_name=folder_name,
        sub="Peptides",
    )

    if cmp_groups is None:
        cmp_groups = data.cmp_groups or [list(data.groups.keys())]

    figures = []

    for _, row in data.psms.iterrows():
        values = []

        for groups in cmp_groups:
            group_vals = pd.Series([
                row[[
                    data.channels[name]
                    for name in data.groups[group]
                    if name in data.channels
                ]]
                for group in groups
            ], index=groups, dtype=object)

            group_vals = pd.Series([
                group[~pd.isnull(group)]
                for group in group_vals
            ], index=group_vals.index, dtype=object)

            group_vals = group_vals[
                group_vals.apply(lambda x: x.shape[0] > 0)
            ]

            # Check normalization group is not null and at least one other
            # group of channels is not null
            if (
                group_vals.shape[0] < 1 or (
                    len(cmp_groups) > 1 and
                    groups[0] not in group_vals.index
                ) or all(
                    group not in group_vals.index
                    for group in groups[1:]
                )
            ):
                continue

            normalize = group_vals.iloc[0].median()

            group_vals = pd.Series([
                group / normalize
                for group in group_vals
            ], index=group_vals.index, dtype=object)

            values.append(group_vals)

        labels = [
            name
            for group in values
            for name in group.index
        ]

        if ax is None:
            _, plot_ax = plt.subplots(
                figsize=(len(labels) * .75, 4),
            )
        else:
            plot_ax = ax

        x = [
            ind
            for ind, l in enumerate(
                j
                for i in values
                for j in i.values
            )
            for k in l
        ]
        y = np.concatenate([
            np.log2(j.astype(float))
            for i in values
            for j in i.values
        ])
        df = pd.DataFrame(
            [
                (
                    np.log2(k),
                    label,
                    "#e19153"
                    if label in set(j[0] for j in cmp_groups) else
                    "#60ae47"
                )
                for i in values
                for label, j in i.iteritems()
                for k in j.values
            ],
            columns=("y", "label", "color"),
        )
        sns.boxplot(
            x="label",
            y="y",
            hue="color",
            data=df,
            ax=plot_ax,
            dodge=False,
            boxprops=dict(alpha=.3),
        )
        sns.swarmplot(
            x=x,
            y=y,
            color=".25",
            ax=plot_ax,
            # size=10,
        )
        plot_ax.axhline(
            0,
            linestyle="--",
            alpha=.25,
        )

        mod_str = row["Modifications"].__str__(prot_index=0)

        plot_ax.set_title(
            title
            if title else
            "{} ({}{})".format(
                row["Sequence"],
                " / ".join(row["Proteins"].genes)[:20],
                (" " + mod_str) if mod_str else "",
            ),
        )
        plot_ax.xaxis.grid(False)

        y_max = y.max()

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

        v = [
            vals
            for group_vals in values
            for vals in group_vals
        ]

        for grp_set in cmp_groups:
            offset = y_max / 10

            for label_a, label_b in itertools.combinations(grp_set, 2):
                if label_a not in labels or label_b not in labels:
                    continue

                index_a = labels.index(label_a)
                index_b = labels.index(label_b)

                values_a = v[index_a]
                values_b = v[index_b]

                if values_a.shape[0] < 2 or values_b.shape[0] < 2:
                    continue

                pval = ttest_ind(
                    values_a.values,
                    values_b.values,
                ).pvalue

                if pval < 0.05:
                    plot_ax.annotate(
                        "",
                        xy=(
                            index_a,
                            y_max + offset,
                        ),
                        xytext=(
                            index_b,
                            y_max + offset,
                        ),
                        xycoords='data',
                        textcoords='data',
                        arrowprops=dict(
                            arrowstyle="-",
                            ec='#000000',
                        ),
                    )
                    plot_ax.text(
                        x=np.mean([index_a, index_b]),
                        y=y_max + offset + y_max / 40,
                        s=stars(pval),
                        horizontalalignment='center',
                        verticalalignment='center',
                    )
                    offset += y_max / 10

            plot_ax.set_ylim(
                bottom=plot_ax.get_ylim()[0],
                top=max([
                    plot_ax.get_ylim()[1],
                    y_max + offset + y_max / 10,
                ]),
            )

        plot_ax.set_xlabel("")
        plot_ax.set_ylabel(
            "{} Signal".format(
                "Relative" if cmp_groups else "Cumulative",
            ),
        )
        plot_ax.get_legend().set_visible(False)

        plot_ax.set_yticklabels(
            ["{:.2f}".format(i) for i in np.power(2, plot_ax.get_yticks())],
        )

        plot_ax.set_xticklabels(
            labels,
            rotation=45,
            horizontalalignment="right",
        )

        plot_ax.get_figure().savefig(
            os.path.join(
                folder_name,
                re.sub(
                    "[?/]",
                    "_",
                    "{} - {} - groups.png".format(
                        "+".join(row["Proteins"].genes)[:50],
                        row["Sequence"],
                    ),
                ),
            ),
            bbox_inches="tight",
            dpi=pyp.DEFAULT_DPI,
            transparent=True,
        )

        figures.append((plot_ax.get_figure(), plot_ax))

    return figures


def plot_all(
    data,
    individual=True,
    between=True,
    cmp_groups=None,
    folder_name=None,
):
    """
    Runs :func:`.plot` and :func:`.plot_group` for all peptides in a data set.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.DataSet`
    figsize : tuple of (int, int), optional
    cmp_groups : list of tuple, optional
    folder_name : str, optional

    Returns
    -------
    figs : list of :class:`matplotlib.figure.Figure`
    """
    folder_name = pyp.utils.make_folder(
        data=data,
        folder_name=folder_name,
        sub="Peptides",
    )

    figures = []

    figures += plot(
        data,
        folder_name=folder_name,
    )

    figures += plot_group(
        data,
        cmp_groups=cmp_groups,
        folder_name=folder_name,
    )

    return figures
