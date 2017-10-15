# -*- coding: UTF-8 -*-
"""
This module provides functionality for data set analysis.

Functions include volcano plots, sorted tables, and plotting sequence levels.
"""

from __future__ import division

# Built-ins
import itertools
import logging
import os
import re

# IPython
# from IPython.display import display

# Core data analysis libraries
from matplotlib import pyplot as plt
# import matplotlib.patches as patches
import matplotlib_venn as mv
import numpy as np
import pandas as pd
import seaborn as sns
# import scipy
from scipy.stats import pearsonr, spearmanr, ttest_ind
# from scipy.stats.mstats import mquantiles
# from scipy.spatial import distance
# from scipy.cluster import hierarchy
# import sklearn
# from sklearn.cluster import KMeans, MiniBatchKMeans

# Misc extras
from adjustText.adjustText import adjust_text
# import fastcluster
# import somoclu
# import uniprot

from . import fetch_data, utils


LOGGER = logging.getLogger("pyproteome.analysis")


def snr_table(
    data,
    sort="p-value",
    folder_name=None, csv_name=None,
):
    """
    Show a table of fold changes and p-values.

    Parameters
    ----------
    data : :class:`DataSet<pyproteome.data_sets.DataSet>`
    sort : str, optional
    folder_name : str, optional
    csv_name : str, optional
    """
    if folder_name is None:
        folder_name = data.name

    if csv_name is None:
        csv_name = "{}-{}.csv".format(
            re.sub("[ ></]", "_", data.name),
            re.sub("[ ></]", "_", data.enrichment),
        )

    utils.make_folder(folder_name)

    csv_name = os.path.join(folder_name, csv_name)

    psms = data.psms[
        [
            "Proteins", "Sequence", "Modifications",
            "Fold Change", "p-value", "Validated",
        ]
    ].copy()
    psms["Sequence"] = [
        "{} ({})".format(row["Sequence"], row["Modifications"])
        for _, row in psms.iterrows()
    ]
    if sort == "Fold Change":
        psms["Fold Change-Sort"] = psms["Fold Change"].apply(
            lambda x: max([x, 1 / x])
        )
        psms.sort_values("Fold Change-Sort", inplace=True, ascending=False)
        psms.drop("Fold Change-Sort", axis=1, inplace=True)
    else:
        psms.sort_values(sort, inplace=True, ascending=True)

    psms.drop("Modifications", axis=1, inplace=True)

    if csv_name:
        psms.to_csv(csv_name)

    # back_colors = {
    #     True: "#BBFFBB",  # light green
    #     False: "#FFBBBB",  # light red
    # }

    if psms.empty:
        return psms

    # return psms.style.apply(  # Color validated rows
    #     lambda row: [
    #         "background-color: " + back_colors[row["Validated"]]
    #         for _ in row
    #     ],
    #     axis=1,
    # )
    return psms.style.set_table_styles(  # Hide index and "Validated" columns
        [
            {"selector": "th:first-child", "props": [("display", "none")]},
            {"selector": "td:last-child", "props": [("display", "none")]},
            {"selector": "th:last-child", "props": [("display", "none")]},
            {"selector": "*", "props": [("text-align", "left")]},
        ]
    )


def _rewrite_enrichments(lst):
    mapping = {
        "MPM2": "pST",
        "MAPK / CDK": "pST",
        "MAPK + CDK": "pST",
        "NTA": "pST",
        "IMAC": "pST",
    }
    return "+".join(
        sorted(
            set(
                mapping.get(i, i)
                for i in lst
            )
        )
    )


def write_full_tables(datas, folder_name="All", out_name="Full Data.xlsx"):
    """
    Write a full list of data sets to a single .xlsx file.

    Parameters
    ----------
    datas : list of :class:`DataSet<pyproteome.data_sets.DataSet>`
    folder_name : str, optional
    out_name : str, optional
    """

    utils.make_folder(folder_name)

    if folder_name is not None:
        out_name = os.path.join(folder_name, out_name)

    writer = pd.ExcelWriter(out_name, engine="xlsxwriter")

    for data in datas:
        channels = [
            (name, data.channels[name])
            for group in data.groups.values()
            for name in group
            if name in data.channels and
            data.channels[name] in data.psms.columns
        ]
        df = data.psms[
            [
                "Proteins", "Sequence",
                "Fold Change", "p-value", "Validated",
            ] + [
                chan
                for _, chan in channels
            ]
        ].copy()

        df.rename(
            columns={
                chan: name
                for name, chan in channels
            },
            inplace=True,
        )
        df.insert(
            2, "Modifications",
            df["Sequence"].apply(
                lambda x: str(x.modifications)
            ),
        )
        df["Sequence"] = df["Sequence"].apply(str)
        df.sort_values("p-value", inplace=True, ascending=True)

        ws_name = re.sub(
            "/",
            "+",
            "{} - {} - {}".format(
                data.name,
                data.tissue,
                _rewrite_enrichments(data.enrichments),
            ),
        )
        df.to_excel(
            writer,
            sheet_name=ws_name,
            index=False,
        )

        ws = writer.sheets[ws_name]
        ws.freeze_panes(1, 0)
        ws.set_column(0, 0, 60)
        ws.set_column(1, 1, 30)
        ws.set_column(2, 2, 20)
        ws.set_column(3, 3, 12)
        ws.set_column(4, 4, 12)

    writer.save()


def _remove_lesser_dups(labels, compress_sym=False):
    new_labels = []

    for index, (p, change, label, color, highlight) in enumerate(labels):
        for o_index, (o_p, o_change, o_label, _, _) in enumerate(labels):
            if index == o_index:
                continue
            if label != o_label:
                continue
            if not compress_sym and (change < 0) != (o_change < 0):
                continue

            mul = 1 if change >= 0 else -1
            o_mul = 1 if o_change >= 0 else -1

            if (
                p + mul * change < o_p + o_mul * o_change
            ) or (
                (
                    p + mul * change <= o_p + o_mul * o_change
                ) and
                index < o_index
            ):
                break
        else:
            new_labels.append((p, change, label, color, highlight))

    return new_labels


def _get_color(txt, x, y):
    return "#BFEE90" if x > 0 else "#FFC1C1"


def volcano_plot(
    data,
    group_a=None,
    group_b=None,
    pval_cutoff=0.05, fold_cutoff=1.25,
    fold_and_p=True,
    xminmax=None, yminmax=None,
    options=None,
    folder_name=None, title=None,
    figsize=(12, 10),
    adjust=True,
    compress_dups=True,
    compress_sym=False,
    show_xlabel=True,
    show_ylabel=True,
    show_title=True,
    full_site_labels=False,
):
    """
    Display a volcano plot of data.

    This plot inclues the fold-changes and p-values associated with said
    changes.

    Parameters
    ----------
    data : :class:`DataSet<pyproteome.data_sets.DataSet>`
    group_a : str or list of str, optional
    group_b : str or list of str, optional
    pval_cutoff : float, optional
    fold_cutoff : float, optional
    options: dict of (str, list), optional
    folder_name : str, optional
    title : str, optional
    figsize : tuple of float, float
    adjust : bool, optional
    compress_dups : bool, optional
    full_site_labels : bool, optional
    """
    (channels_a, channels_b), (label_a, label_b), _ = data.get_groups(
        group_a=group_a,
        group_b=group_b,
    )

    if group_a and group_b:
        data.update_group_changes(group_a=group_a, group_b=group_b)

    if not folder_name:
        folder_name = data.name

    options = options or {}

    highlight = options.get('highlight', {})
    hide = options.get('hide', {})
    show = options.get('show', {})
    edgecolors = options.get('edgecolors', {})
    rename = options.get('rename', {})

    log_pval_cutoff = -np.log10(pval_cutoff)

    utils.make_folder(folder_name)

    if not title:
        title = "{} - {}".format(
            data.tissue,
            data.enrichment,
        )

    if title:
        file_name = re.sub("[ ></?]", "_", title) + "_Volcano.png"

        if folder_name:
            file_name = os.path.join(folder_name, file_name)

    upper_fold = np.log2(fold_cutoff)
    lower_fold = -upper_fold

    # Calculate the Fold-Change / p-values
    pvals = []
    changes = []
    colors = []
    labels = []

    for _, row in data.psms.dropna(
        subset=["p-value", "Fold Change"],
    ).iterrows():
        color = "grey"
        row_pval = -np.log10(row["p-value"])
        row_change = np.log2(row["Fold Change"])

        row_label = " / ".join(sorted(row["Proteins"].genes))
        old_row_label, old_re_row_label = None, None

        if (
            full_site_labels and
            len(list(row["Modifications"].skip_labels_iter())) > 0
        ):
            old_row_label = row_label
            old_re_row_label = rename.get(old_row_label, old_row_label)
            row_label = " / ".join(
                sorted(
                    "{} {}".format(
                        gene,
                        row["Modifications"].__str__(prot_index=index),
                    )
                    for index, gene in enumerate(row["Proteins"].genes)
                )
            )

        re_row_label = rename.get(row_label, row_label)
        edgecolor = edgecolors.get(
            re_row_label, edgecolors.get(
                row_label,
                edgecolors.get(old_row_label, None),
            )
        )

        if edgecolor is None:
            gene_colors = [
                edgecolors.get(gene, None)
                for gene in row["Proteins"].genes
            ]
            if len(set(gene_colors)) == 1:
                edgecolor = gene_colors[0]

        if (
            np.isnan(row_pval) or
            np.isnan(row_change) or
            np.isinf(row_pval) or
            np.isinf(row_change)
        ):
            continue

        pvals.append(row_pval)
        changes.append(row_change)

        names = [
            row_label, re_row_label, old_row_label, old_re_row_label,
        ]

        if (
            any(
                i in show
                for i in names
            ) or ((
                row_pval > log_pval_cutoff and
                (row_change > upper_fold or row_change < lower_fold)
            ) if fold_and_p else (
                row_pval > log_pval_cutoff or
                (row_change > upper_fold or row_change < lower_fold)
            ))
        ):
            color = "blue"
            highlight_label = any(i in highlight for i in names)

            if (
                not any(i.strip() in hide for i in row_label.split("/")) and
                re_row_label not in hide
            ):
                labels.append((
                    row_pval, row_change,
                    re_row_label, edgecolor, highlight_label,
                ))

        colors.append(color)

    if compress_dups:
        labels = _remove_lesser_dups(labels, compress_sym=compress_sym)

    # Draw the figure
    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(changes, pvals, c=colors)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if xminmax:
        ax.set_xlim(xmin=xminmax[0], xmax=xminmax[1])
    else:
        ax.set_xlim(
            xmin=np.floor(min(changes) * 2) / 2,
            xmax=np.ceil(max(changes) * 2) / 2,
        )

    if yminmax:
        ax.set_ylim(bottom=yminmax[0], top=yminmax[1])
    else:
        ax.set_ylim(bottom=-0.1)

    ax.set_xticks(
        list(
            sorted(
                tick
                for tick in ax.get_xticks()
                if tick < lower_fold or tick > upper_fold
            ) + [lower_fold, upper_fold]
        )
    )
    ax.set_yticks(
        list(
            sorted(
                [
                    tick
                    for tick in ax.get_yticks()
                    if "{}".format(np.power(1/10, tick)).strip("0.")[:1] in
                    ["1", "5"] and
                    tick > log_pval_cutoff
                ] + [log_pval_cutoff]
            )
        )
    )
    ax.set_xticklabels(
        [
            "{:.3}".format(i)
            for i in np.exp2(ax.get_xticks())
        ],
        fontsize=20,
    )
    ax.set_yticklabels(
        [
            ("{:.3}" if i > 0.005 else "{:.0e}").format(i)
            for i in np.power(1/10, ax.get_yticks())
        ],
        fontsize=20,
    )

    if show_xlabel:
        ax.set_xlabel(
            "Fold Change {} (n={}) / {} (n={})".format(
                label_a,
                len(channels_a),
                label_b,
                len(channels_b),
            ),
            fontsize=20,
        )

    if show_ylabel:
        ax.set_ylabel(
            "p-value",
            fontsize=20,
        )

    if not np.isnan(log_pval_cutoff):
        ax.axhline(
            log_pval_cutoff,
            color="r", linestyle="dashed", linewidth=0.5,
        )

    if abs(fold_cutoff - 1) > 0.01:
        ax.axvline(upper_fold, color="r", linestyle="dashed", linewidth=0.5)
        ax.axvline(lower_fold, color="r", linestyle="dashed", linewidth=0.5)

    # Position the labels
    texts = []
    txt_lim = 100 if full_site_labels else 20

    for y, x, txt, edgecolor, highlight_label in labels:
        text = ax.text(
            x, y,
            txt[:txt_lim] + ("..." if len(txt) > txt_lim else ""),
        )

        if highlight_label:
            text.set_fontsize(20)

        text.set_bbox(
            dict(
                facecolor=_get_color(txt, x, y),
                alpha=1,
                linewidth=0.5 if not edgecolor else 3,
                edgecolor=edgecolor or "black",
                boxstyle="round",
            )
        )

        texts.append(text)

    if adjust:
        adjust_text(
            x=[i[0] for i in labels],
            y=[i[1] for i in labels],
            texts=texts,
            ax=ax,
            lim=500,
            force_text=0.5,
            force_points=0.01,
            arrowprops=dict(arrowstyle="->", relpos=(0, 0), lw=1),
            only_move={
                "points": "y",
                "text": "xy",
            }
        )

    if title:
        if show_title:
            ax.set_title(
                title,
                fontsize=20,
            )

        fig.savefig(
            file_name,
            bbox_inches="tight", dpi=300,
            transparent=True,
        )
        fig.savefig(
            os.path.splitext(file_name)[0] + ".svg",
            bbox_inches="tight", dpi=300,
            transparent=True,
        )

    return fig, ax


def hierarchical_heatmap(
    data,
    baseline_channels=None,
    metric="euclidean", method="centroid",
    row_cluster=True, col_cluster=True
):
    """
    Plot a hierarhically-clustered heatmap of a data set.

    Parameters
    ----------
    data : :class:`DataSet<pyproteome.data_sets.DataSet>`
    baseline_channels : list of str, optional
        List of channels to average and use as baseline for each row.
    metric : str, optional
        Hierarchical clustering distance metric.
    method : str, optional
        Hierarchical clustering method.
    row_cluster : bool, optional
    col_cluster : bool, optional
    """
    channels = [
        channel
        for channels in data.groups.values()
        for channel in channels
    ]
    channel_names = [
        data.channels[i]
        for i in channels
    ]

    if baseline_channels is None:
        baseline_channels = list(data.groups.values())[-1]

    psms = data.psms[
        data.psms["Fold Change"].apply(lambda x: max([x, 1/x])) > 1.5
    ]

    raw = psms[channels].as_matrix()
    raw_baseline = psms[baseline_channels].as_matrix().mean(axis=1)
    raw = np.log2((raw.T / raw_baseline).T)

    sns.clustermap(
        raw,
        method=method,
        metric=metric,
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        xticklabels=channel_names,
        yticklabels=False,
    )


def venn2(data_a, data_b, folder_name=None, filename=None):
    """
    Display a three-way venn diagram between data set sequences.

    Parameters
    ----------
    data_a : :class:`DataSet<pyproteome.data_sets.DataSet>`
    data_b : :class:`DataSet<pyproteome.data_sets.DataSet>`
    data_c : :class:`DataSet<pyproteome.data_sets.DataSet>`
    folder_name : str, optional
    filename : str, optional
    """
    utils.make_folder(folder_name)

    if folder_name and filename:
        filename = os.path.join(folder_name, filename)

    group_a = set(data_a["Sequence"])
    group_b = set(data_b["Sequence"])

    f = plt.figure(figsize=(12, 12))
    v = mv.venn2(
        subsets=(
            len(group_a.difference(group_b)),
            len(group_b.difference(group_a)),
            len(group_a.intersection(group_b)),
        ),
        set_labels=(data_a.tissue, data_b.tissue),
    )

    for label in v.set_labels:
        if label:
            label.set_fontsize(32)

    for label in v.subset_labels:
        if label:
            label.set_fontsize(20)

    if filename:
        f.savefig(filename, transparent=True, dpi=500)


def venn3(data_a, data_b, data_c, folder_name=None, filename=None):
    """
    Display a three-way venn diagram between data set sequences.

    Parameters
    ----------
    data_a : :class:`DataSet<pyproteome.data_sets.DataSet>`
    data_b : :class:`DataSet<pyproteome.data_sets.DataSet>`
    data_c : :class:`DataSet<pyproteome.data_sets.DataSet>`
    folder_name : str, optional
    filename : str, optional
    """
    utils.make_folder(folder_name)

    if folder_name and filename:
        filename = os.path.join(folder_name, filename)

    group_a = set(data_a["Sequence"])
    group_b = set(data_b["Sequence"])
    group_c = set(data_c["Sequence"])

    f = plt.figure(figsize=(12, 12))
    v = mv.venn3(
        subsets=(
            len(group_a.difference(group_b).difference(group_c)),
            len(group_b.difference(group_a).difference(group_c)),
            len(group_a.intersection(group_b).difference(group_c)),
            len(group_c.difference(group_a).difference(group_b)),
            len(group_a.intersection(group_c).difference(group_b)),
            len(group_b.intersection(group_c).difference(group_a)),
            len(group_a.intersection(group_b).intersection(group_c)),
        ),
        set_labels=(data_a.tissue, data_b.tissue, data_c.tissue),
    )

    for label in v.set_labels:
        if label:
            label.set_fontsize(32)

    for label in v.subset_labels:
        if label:
            label.set_fontsize(20)

    if filename:
        f.savefig(filename, transparent=True, dpi=500)


def write_lists(
    data,
    folder_name=None, sorted_name="sorted_list.txt",
    hits_name="hits_list.txt", background_name="back_list.txt",
):
    """
    Write a list of peptides to files.

    Includes peptides sorted by fold change, significantly changing peptides,
    and background peptides with little change.

    Parameters
    ----------
    data : :class:`DataSet<pyproteome.data_sets.DataSet>`
    folder_name : str, optional
    sorted_name : str, optional
    hits_name : str, optional
    background_name : str, optional
    """
    if folder_name is None:
        folder_name = data.name

    utils.make_folder(folder_name)

    if folder_name:
        sorted_name = os.path.join(folder_name, sorted_name)
        hits_name = os.path.join(folder_name, hits_name)
        background_name = os.path.join(folder_name, background_name)

    change_psms = data.psms.copy()
    change_psms["Fold Change"] = np.maximum.reduce(
        [
            change_psms["Fold Change"],
            1 / change_psms["Fold Change"],
        ]
    )

    with open(sorted_name, "w") as f:
        f.write(
            "\n".join(
                i.accessions[0]
                for i in change_psms.sort(
                    "Fold Change",
                    ascending=False,
                )["Proteins"].drop_duplicates(keep="first")
            )
        )
    with open(hits_name, "w") as f:
        f.write(
            "\n".join(
                i.accessions[0]
                for i in data.filter(
                    fold_cutoff=1.3,
                ).psms["Proteins"].drop_duplicates(keep="first")
            )
        )

    with open(background_name, "w") as f:
        f.write(
            "\n".join(
                i.accessions[0]
                for i in data.psms["Proteins"].drop_duplicates(keep="first")
            )
        )


def plot_sequence_between(
    data, sequences, cmp_groups=None,
):
    """
    Plot the levels of a sequence across each group.

    Parameters
    ----------
    data : :class:`DataSet<pyproteome.data_sets.DataSet>`
    sequences : list of str
    """
    channels = [
        [
            data.channels[channel_name]
            for channel_name in group
            if channel_name in data.channels
        ]
        for group in data.groups.values()
    ]

    psms = data.psms.copy()
    psms["Seq Str"] = psms["Sequence"].apply(str)
    psms = psms[psms["Seq Str"].isin(sequences)]

    values = [
        psms[channel].as_matrix().ravel()
        for channel in channels
    ]
    values = [
        row[~np.isnan(row)]
        for row in values
    ]
    values, labels = zip(*[
        (value, label)
        for label, value in zip(data.groups.keys(), values)
        if value.any() and
        (cmp_groups is None or any(label in i for i in cmp_groups))
    ])

    means = np.array([row.mean() for row in values])
    errs = np.array([row.std() for row in values])

    if cmp_groups:
        div = np.array([
            means[[
                labels.index(group[0])
                for group in cmp_groups
                if label in group
            ][0]]
            for label in labels
        ])

        means /= div
        errs /= div

    f, ax = plt.subplots()

    indices = np.arange(len(means))
    bar_width = .35
    ax.bar(
        bar_width * 1.5 + indices,
        means,
        bar_width,
        yerr=errs,
        error_kw=dict(ecolor='k', lw=2, capsize=10, capthick=2),
    )

    ax.set_ylabel(
        "{} Signal".format(
            "Relative" if cmp_groups else "Cumulative",
        ),
        fontsize=20,
    )
    ax.ticklabel_format(style="plain")

    for label in ax.get_yticklabels():
        label.set_fontsize(14)

    ax.set_xticks(indices + bar_width * 1.5)
    ax.set_xticklabels(labels, fontsize=16)

    title = "{} - {}".format(
        " / ".join(sequences),
        " / ".join(sorted(
            gene
            for row in psms["Proteins"]
            for gene in row.genes
        ))
    )
    ax.set_title(title, fontsize=20)
    ax.xaxis.grid(False)

    if not means.any():
        return f

    y_max = np.max(means + errs)

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
                    index_a + bar_width * 1.5,
                    y_max * (1.05 + offset * 0.1)
                ),
                xytext=(
                    index_b + bar_width * 1.5,
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
                x=np.mean([index_a, index_b]) + bar_width * 1.5,
                y=y_max * (1.07 + offset * 0.1),
                s=stars(pval),
                horizontalalignment='center',
                verticalalignment='center',
            )

            offset += 1

    return f


def plot_sequence(
    data, sequence,
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
    channels = [
        data.channels[channel_name]
        for channel_name in channel_names
    ]

    psms = data.psms[data.psms["Sequence"] == sequence]

    values = psms[channels].as_matrix()

    mask = ~np.isnan(values).all(axis=0)
    channel_names = list(np.array(channel_names)[mask])
    values = values[:, mask]

    f, ax = plt.subplots(figsize=figsize)

    for i in range(values.shape[0]):
        indices = np.arange(len(values[i]))
        bar_width = .35
        ax.bar(indices, values[i], bar_width)
        ax.set_xticks(indices)
        ax.set_xticklabels(
            channel_names,
            fontsize=20,
            rotation=45,
            horizontalalignment="right",
        )

    ax.set_title(sequence + (" - {}".format(title) if title else ""))

    ax.set_ylabel(
        "Cummulative Intensity" +
        (" (Normalized)" if data.intra_normalized else "")
    )

    ax.title.set_fontsize(28)
    ax.yaxis.label.set_fontsize(20)

    return f


def correlate_data_sets(
    data1, data2, folder_name=None, filename=None,
    adjust=True, label_cutoff=1.5,
):
    """
    Plot the correlation between peptides levels in two different data sets.

    Parameters
    ----------
    data1 : :class:`DataSet<pyproteome.data_sets.DataSet>`
    data2 : :class:`DataSet<pyproteome.data_sets.DataSet>`
    folder_name : str, optional
    filename : str, optional
    """
    if not folder_name:
        folder_name = "All" if data1.name != data2.name else data1.name

    utils.make_folder(folder_name)

    if folder_name and filename:
        filename = os.path.join(folder_name, filename)

    merged = pd.merge(
        data1.psms, data2.psms,
        on="Sequence",
    ).dropna()

    f, ax = plt.subplots()
    ax.scatter(
        np.log2(merged["Fold Change_x"]),
        np.log2(merged["Fold Change_y"]),
    )

    label_cutoff = np.log2(label_cutoff)
    texts, x_s, y_s = [], [], []

    for index, row in merged.iterrows():
        x = row["Fold Change_x"]
        y = row["Fold Change_y"]
        ratio = np.log2(x / y)
        x, y = np.log2(x), np.log2(y)

        if ratio < label_cutoff and ratio > - label_cutoff:
            continue

        x_s.append(x)
        y_s.append(y)

        txt = " / ".join(row["Proteins_x"].genes)
        txt = txt[:20] + ("..." if len(txt) > 20 else "")

        text = ax.text(
            x, y, txt,
        )

        text.set_bbox(
            dict(
                color="lightgreen" if ratio < 0 else "pink",
                alpha=0.8,
            )
        )
        texts.append(text)

    if adjust:
        adjust_text(
            x=x_s,
            y=y_s,
            texts=texts,
            ax=ax,
            lim=400,
            force_text=0.1,
            force_points=0.1,
            arrowprops=dict(arrowstyle="->", relpos=(0, 0), lw=1),
            only_move={
                "points": "y",
                "text": "xy",
            }
        )
    else:
        for j, text in enumerate(texts):
            a = ax.annotate(
                text.get_text(),
                xy=text.get_position(),
                xytext=text.get_position(),
                arrowprops=dict(arrowstyle="->", relpos=(0, 0), lw=1),
            )
            a.__dict__.update(text.__dict__)
            a.draggable()
            texts[j].remove()

    min_x = min(np.log2(merged["Fold Change_x"]))
    max_x = max(np.log2(merged["Fold Change_x"]))
    ax.plot([min_x, max_x], [min_x, max_x], "--")
    ax.plot(
        [min_x, max_x],
        [min_x + label_cutoff, max_x + label_cutoff],
        color="lightgreen",
        linestyle=":",
    )
    ax.plot(
        [min_x, max_x],
        [min_x - label_cutoff, max_x - label_cutoff],
        color="pink",
        linestyle=":",
    )

    name1 = "{} - {}".format(data1.name, ", ".join(sorted(data1.tissues)))
    name2 = "{} - {}".format(data2.name, ", ".join(sorted(data2.tissues)))

    ax.set_xlabel("$log_2$ Fold Change -- {}".format(name1))
    ax.set_ylabel("$log_2$ Fold Change -- {}".format(name2))

    pear_corr = pearsonr(merged["Fold Change_x"], merged["Fold Change_y"])
    spear_corr = spearmanr(merged["Fold Change_x"], merged["Fold Change_y"])

    ax.set_title(
        (
            r"Pearson's: $\rho$={:.2f}, "
            r"Spearman's: $\rho$={:.2f}"
        ).format(pear_corr[0], spear_corr[0])
    )

    if filename:
        f.savefig(filename, transparent=True, dpi=500)


def find_tfs(data, folder_name=None, csv_name=None):
    """
    Scan over a data set to find proteins annotated as transcription factors.

    Parameters
    ----------
    data : :class:`DataSet<pyproteome.data_sets.DataSet>`
    folder_name : str, optional
    csv_name : str, optional
    """
    if folder_name is None:
        folder_name = data.name

    utils.make_folder(folder_name)

    if csv_name is None:
        csv_name = "Changing TFs.csv"

    if folder_name and csv_name:
        csv_name = os.path.join(folder_name, csv_name)

    def _is_tf(prots):
        go_terms = (
            "DNA binding",
            "double-stranded DNA binding",
            "transcription factor binding",
            "transcription, DNA-templated",
        )
        return any(
            go_term in go
            for prot in prots
            for go in fetch_data.get_uniprot_data(prot.accession).get(
                "go", [],
            )
            for go_term in go_terms
        )

    tfs = data.psms[data.psms["Proteins"].apply(_is_tf)].copy()
    tfs.sort_values(by="Fold Change", ascending=False, inplace=True)

    if csv_name:
        tfs[["Proteins", "Sequence", "Modifications", "Fold Change"]].to_csv(
            csv_name,
            index=False,
        )

    tfs = tfs[["Proteins", "Sequence", "Fold Change", "p-value"]]

    return tfs.style.set_table_styles(  # Hide index and "Validated" columns
        [
            {"selector": "th:first-child", "props": [("display", "none")]},
            {"selector": "*", "props": [("text-align", "left")]},
        ]
    )


def spearmanr_nan(a, b, min_length=5):
    mask = ~np.array(
        [np.isnan(i) or np.isnan(j) for i, j in zip(a, b)],
        dtype=bool,
    )

    if mask.sum() < min_length:
        return spearmanr(np.nan, np.nan)

    return spearmanr(a[mask], b[mask])


def correlate_signal(
    data, signal,
    pval_cutoff=0.05, fold_cutoff=1.25,
    options=None,
    folder_name=None, title=None,
    scatter_colors=None,
    scatter_symbols=None,
    figsize=(12, 10),
    xlabel="",
):

    options = options or {}

    highlight = options.get('highlight', {})
    hide = options.get('hide', {})
    edgecolors = options.get('edgecolors', {})
    rename = options.get('rename', {})
    scatter_colors = options.get('scatter_colors', {})
    scatter_symbols = options.get('scatter_symbols', {})

    cp = data.copy()

    signal_groups = [
        label
        for label, group in cp.groups.items()
        for chan in group
        if chan in cp.channels.keys() and
        chan in signal.columns
    ]

    signal_chans = [
        chan
        for group in cp.groups.values()
        for chan in group
        if chan in cp.channels.keys() and
        chan in signal.columns
    ]
    data_chans = [
        data.channels[chan]
        for chan in signal_chans
    ]

    corr = [
        spearmanr_nan(
            row[data_chans].as_matrix().ravel(),
            signal[signal_chans].as_matrix().ravel(),
        )
        for _, row in cp.psms.iterrows()
    ]

    cp.psms["Correlation"] = [i.correlation for i in corr]
    cp.psms["corr p-value"] = [i.pvalue for i in corr]

    f_corr, ax = plt.subplots(figsize=figsize)
    x, y, colors = [], [], []
    sig_x, sig_y, sig_labels = [], [], []

    for _, row in cp.psms.iterrows():
        if (
            row["Correlation"] == 0 or
            row["corr p-value"] == 0 or
            np.isinf(row["Correlation"]) or
            np.isnan(row["Correlation"]) or
            np.isinf(-np.log10(row["corr p-value"])) or
            np.isnan(-np.log10(row["corr p-value"]))
        ):
            continue

        x.append(row["Correlation"])
        y.append(-np.log10(row["corr p-value"]))

        if row["corr p-value"] < pval_cutoff:
            sig_x.append(row["Correlation"])
            sig_y.append(-np.log10(row["corr p-value"]))
            sig_labels.append(" / ".join(sorted(row["Proteins"].genes)))

        colors.append(
            "blue"
            if row["corr p-value"] < pval_cutoff else
            "{:.2f}".format(
                max([len(row[data_chans].dropna()) / len(data_chans) - .25, 0])
            )
        )

    ax.scatter(x, y, c=colors)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    texts = []
    for xs, ys, txt in zip(sig_x, sig_y, sig_labels):
        if any(i.strip() in hide for i in txt.split("/")):
            continue

        txt = rename.get(txt, txt)

        text = ax.text(
            xs, ys,
            txt[:20] + ("..." if len(txt) > 20 else ""),
        )

        if txt in highlight:
            text.set_fontsize(20)

        text.set_bbox(
            dict(
                facecolor="lightgreen" if xs > 0 else "pink",
                alpha=1,
                linewidth=0.5 if txt not in edgecolors else 3,
                edgecolor=edgecolors.get(txt, "black"),
                boxstyle="round",
            )
        )

        texts.append(text)

    adjust_text(
        x=sig_x,
        y=sig_y,
        texts=texts,
        ax=ax,
        lim=400,
        force_text=0.3,
        force_points=0.01,
        arrowprops=dict(arrowstyle="->", relpos=(0, 0), lw=1),
        only_move={
            "points": "y",
            "text": "xy",
        }
    )

    ax.set_xlabel(
        "Correlation",
        fontsize=20,
    )
    ax.set_yticklabels(
        "{:.3}".format(i)
        for i in np.power(1/10, ax.get_yticks())
    )
    ax.set_ylabel(
        "p-value",
        fontsize=20,
    )

    cp.psms = cp.psms[cp.psms["corr p-value"] < pval_cutoff]

    f_scatter, axes = plt.subplots(
        int(np.ceil(cp.psms.shape[0] / 3)), 3,
        figsize=(18, cp.psms.shape[0] * 2),
    )

    for index, (ax, (_, row)) in enumerate(
        zip(
            axes.ravel(),
            cp.psms.sort_values("Correlation").iterrows(),
        )
    ):
        for data_chan, sig_chan, sig_group in zip(
            data_chans, signal_chans, signal_groups,
        ):
            ax.scatter(
                x=signal[sig_chan],
                y=row[data_chan],
                facecolors=scatter_colors.get(sig_group, "black"),
                edgecolors="black",
                marker=scatter_symbols.get(sig_group, "o"),
                s=200,
            )

        row_title = " / ".join(str(i.gene) for i in row["Proteins"])

        ax.set_title(
            row_title[:20] + ("..." if len(row_title) > 20 else ""),
            fontsize=28,
            fontweight="bold",
        )

        ax.set_xlabel(
            "{}\n$\\rho = {:.2f}; p = {:.2E}$".format(
                xlabel,
                row["Correlation"],
                row["corr p-value"],
            ),
            fontsize=22,
        )

        row_seq = str(row["Sequence"])
        row_mods = str(row["Modifications"].get_mods([(None, "Phospho")]))

        ax.set_ylabel(
            row_seq[:20] + (
                "..." if len(row_seq) > 20 else ""
            ) + (
                "\n({})".format(
                    row_mods[:20] + ("..." if len(row_mods) > 20 else "")
                ) if row_mods else ""
            ),
            fontsize=20,
        )
        for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(20)

    f_scatter.tight_layout(pad=2)

    return f_corr, f_scatter


def plot_all(
    datas, seqs=None, protein=None, figsize=(16, 8),
    individual=True, between=False, cmp_groups=None,
):
    assert seqs is not None or protein is not None

    if protein:
        seqs = list(
            set(
                str(seq)
                for data in datas
                for seq in data[data["Proteins"] == protein]["Sequence"]
            )
        )

    if isinstance(seqs, str):
        seqs = [seqs]

    for seq in seqs:
        for data in datas:
            prot = " / ".join(
                gene
                for i in data[data["Sequence"] == seq]["Proteins"]
                for gene in i.genes
            )
            if individual:
                f = plot_sequence(
                    data, seq,
                    title="{}-{}".format(
                        prot,
                        data.name,
                        ":".join(data.tissues),
                    ),
                    figsize=figsize,
                )
                f.savefig(
                    re.sub(
                        "[?/]",
                        "_",
                        "{}-{}-{}.png".format(
                            prot,
                            ",".join(data.tissues),
                            data.name,
                        ),
                    ),
                    bbox_inches="tight", dpi=300,
                    transparent=True,
                )

            if between:
                f = plot_sequence_between(data, [seq], cmp_groups=cmp_groups)
                f.savefig(
                    re.sub(
                        "[?/]",
                        "_",
                        "{}-{}-{}-between.png".format(
                            prot,
                            ",".join(data.tissues),
                            data.name,
                        ),
                    ),
                    bbox_inches="tight", dpi=300,
                    transparent=True,
                )


def plot_volcano_filtered(data, f, **kwargs):
    data = data.copy()
    data.update_group_changes(
        group_a=kwargs.get("group_a", None),
        group_b=kwargs.get("group_b", None),
    )

    d = data.filter(**f)

    changes = []
    pvals = []

    for _, row in data.psms.iterrows():
        row_pval = -np.log10(row["p-value"])
        row_change = np.log2(row["Fold Change"])

        if (
            np.isnan(row_pval) or
            np.isnan(row_change) or
            np.isinf(row_pval) or
            np.isinf(row_change)
        ):
            continue

        pvals.append(row_pval)
        changes.append(row_change)

    xminmax = kwargs.pop(
        "xminmax",
        (
            np.floor(min(changes) * 2) / 2,
            np.ceil(max(changes) * 2) / 2,
        ),
    )
    yminmax = kwargs.pop(
        "yminmax",
        (-0.1, np.ceil(max(pvals))),
    )

    f, ax = volcano_plot(
        d,
        xminmax=xminmax,
        yminmax=yminmax,
        **kwargs
    )

    changes, pvals = zip(*[
        (x, y)
        for x, y in zip(changes, pvals)
        if x > xminmax[0] and x < xminmax[1] and
        y > yminmax[0] and y < yminmax[1]
    ])

    ax.scatter(changes, pvals, c="lightblue", zorder=0, alpha=0.3)

    f.savefig(
        os.path.join(
            data.name,
            re.sub(r"[ ></\?]", "_", kwargs.get("title", "Filtered")) +
            "_Volcano.png",
        ),
        bbox_inches="tight", dpi=300,
        transparent=True,
    )
    return f, ax
