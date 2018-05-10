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
    data, f=None,
    title=None,
    folder_name=None,
    figsize=None,
):
    """
    Plot the levels of a sequence across multiple channels.

    Parameters
    ----------
    data : :class:`DataSet<pyproteome.data_sets.DataSet>`
    sequence : str or :class:`Sequence<pyproteome.sequence.Sequence>`
    title : str, optional
    figsize : tuple of (int, int), optional
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

    if figsize is None:
        figsize = (len(channels) / 2, 6 / 2)

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

        sns.barplot(
            x="name",
            y="val",
            hue="group",
            data=df,
            ax=ax,
            dodge=False,
        )
        ax.set_xticklabels(
            ax.get_xticklabels(),
            # fontsize=20,
            rotation=45,
            horizontalalignment="right",
        )
        ax.set_xlabel("")
        ax.get_legend().set_title("")

        mod_str = row["Modifications"].__str__(prot_index=0)

        ax.set_title(
            title
            if title else
            "{} ({}{})".format(
                seq,
                " / ".join(row["Proteins"].genes)[:20],
                (" " + mod_str) if mod_str else "",
            ),
            fontsize=28,
        )

        if data.cmp_groups:
            ylabel = "Relative Intensity"
        else:
            ylabel = (
                "Cummulative Intensity" +
                (" (Normalized)" if data.intra_normalized else "")
            )

        ax.set_ylabel(
            ylabel,
            fontsize=20,
        )

        fig.savefig(
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

        figures.append((fig, ax))

    return figures


def plot_group(
    data,
    f=None,
    cmp_groups=None,
    title=None,
    folder_name=None,
    figsize=(8, 4),
):
    """
    Plot the levels of a sequence across each group.

    Parameters
    ----------
    data : :class:`DataSet<pyproteome.data_sets.DataSet>`
    sequences : list of str
    """
    folder_name = pyp.utils.make_folder(
        data=data,
        folder_name=folder_name,
        sub="Peptides",
    )

    if cmp_groups is None:
        cmp_groups = data.cmp_groups or [list(data.groups.keys())]

    if f:
        data = data.filter(f)

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

        # print(values)
        fig, ax = plt.subplots(figsize=figsize)

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
        sns.boxplot(
            x=x,
            y=y,
            ax=ax,
            boxprops=dict(alpha=.3),
        )
        sns.swarmplot(
            x=x,
            y=y,
            color=".25",
            ax=ax,
            size=10,
        )
        ax.axhline(
            0,
            linestyle="--",
            alpha=.25,
        )

        means = [
            i.mean()
            for vals in values
            for i in vals
        ]
        means = [
            np.log2(i.mean())
            for vals in values
            for i in vals
        ]

        errs = [
            np.log2(i.std()) if len(i) > 1 else 0
            for vals in values
            for i in vals
        ]

        ax.set_ylabel(
            "{} Signal".format(
                "Relative" if cmp_groups else "Cumulative",
            ),
            fontsize=20,
        )

        ax.set_yticklabels(
            ["{:.2f}".format(i) for i in np.power(2, ax.get_yticks())],
            fontsize=20,
        )

        ax.set_xticklabels(
            labels,
            fontsize=20,
            rotation=45,
            horizontalalignment="right",
        )

        mod_str = row["Modifications"].__str__(prot_index=0)

        ax.set_title(
            title
            if title else
            "{} ({}{})".format(
                row["Sequence"],
                " / ".join(row["Proteins"].genes)[:20],
                (" " + mod_str) if mod_str else "",
            ),
            fontsize=32,
        )
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

        v = [
            vals
            for group_vals in values
            for vals in group_vals
        ]

        for grp_set in cmp_groups:
            offset = 0

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
                    values_a.as_matrix(),
                    values_b.as_matrix(),
                ).pvalue

                # if pval < 0.05:
                #     ax.set_ylim(
                #         ymin=ax.get_ylim()[0],
                #         ymax=max([
                #             ax.get_ylim()[1],
                #             y_max * (1 + offset * 0.1 + 0.15),
                #         ]),
                #     )
                #     ax.annotate(
                #         "",
                #         xy=(
                #             index_a,
                #             y_max * (1.05 + offset * 0.1)
                #         ),
                #         xytext=(
                #             index_b,
                #             y_max * (1.05 + offset * 0.1)
                #         ),
                #         xycoords='data',
                #         textcoords='data',
                #         arrowprops=dict(
                #             arrowstyle="-",
                #             ec='#000000',
                #         ),
                #     )
                #     ax.text(
                #         x=np.mean([index_a, index_b]),
                #         y=y_max * (1.07 + offset * 0.1),
                #         s=stars(pval),
                #         horizontalalignment='center',
                #         verticalalignment='center',
                #     )
                #
                #     offset += 1

        fig.savefig(
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

        figures.append((fig, ax))

    return figures


def plot_all(
    data,
    f=None,
    figsize=(16, 8),
    individual=True,
    between=True,
    cmp_groups=None,
    folder_name=None,
):
    folder_name = pyp.utils.make_folder(
        data=data,
        folder_name=folder_name,
        sub="Peptides",
    )

    figures = []

    if individual:
        figures += plot(
            data,
            f=f,
            figsize=figsize,
            folder_name=folder_name,
        )

    if between:
        figures += plot_group(
            data,
            f=f,
            cmp_groups=cmp_groups,
            folder_name=folder_name,
        )

    return figures


def plot_together(
    data,
    folder_name=None,
    only=True,
    **kwargs
):
    folder_name = pyp.utils.make_folder(
        data=data,
        folder_name=folder_name,
        sub="Peptides",
    )

    figures = plot_group(
        data,
        folder_name=folder_name,
        **kwargs
    )

    return figures
