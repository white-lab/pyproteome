"""
This module provides functionality for data set analysis.

Functions include volcano plots, sorted tables, and plotting sequence levels.
"""

from __future__ import division

# Built-ins
import logging
import os
import re

# IPython
from IPython.display import display

# Core data analysis libraries
from matplotlib import pyplot as plt
# import matplotlib.patches as patches
import matplotlib_venn as mv
import numpy as np
import pandas as pd
import seaborn as sns
# import scipy
from scipy.stats import ttest_ind, pearsonr, spearmanr
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
    Show a signal to noise table.

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
    ]
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

    back_colors = {
        True: "#BBFFBB",  # light green
        False: "#FFBBBB",  # light red
    }

    return psms.style.apply(  # Color validated rows
        lambda row: [
            "background-color: " + back_colors[row["Validated"]]
            for _ in row
        ],
        axis=1,
    ).set_table_styles(  # Hide index and "Validated" columns
        [
            {"selector": "th:first-child", "props": [("display", "none")]},
            {"selector": "td:last-child", "props": [("display", "none")]},
            {"selector": "th:last-child", "props": [("display", "none")]},
        ]
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
        df = data.psms[
            [
                "Proteins", "Sequence",
                "Fold Change", "p-value", "Validated",
            ]
        ]
        df.insert(
            2, "Modifications",
            df["Sequence"].apply(
                lambda x: str(x.modifications)
            ),
        )
        df["Sequence"] = df["Sequence"].apply(str)
        df.sort_values("p-value", inplace=True, ascending=True)

        ws_name = "{}-{}".format(data.name, data.enrichment)
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


def _remove_lesser_dups(pvals, changes, labels):
    new_pvals, new_changes, new_labels = [], [], []

    for index, (p, change, label) in enumerate(zip(pvals, changes, labels)):
        for o_index, o_label in enumerate(labels):
            if index == o_index:
                continue
            if (
                label == o_label and
                p + change < pvals[o_index] + changes[o_index]
            ):
                break
        else:
            new_pvals.append(p)
            new_changes.append(change)
            new_labels.append(label)

    return new_pvals, new_changes, new_labels


def volcano_plot(
    data,
    pval_cutoff=1.3, fold_cutoff=1.2,
    highlight=None,
    folder_name=None, title=None,
    figsize=(12, 10),
    adjust=True,
    compress_dups=True,
):
    """
    Display a volcano plot of data.

    This plot inclues the fold-changes and p-values associated with said
    changes.

    Parameters
    ----------
    data : :class:`DataSet<pyproteome.data_sets.DataSet>`
    pval_cutoff : float, optional
    fold_cutoff : float, optional
    highlight : list, optional
    folder_name : str, optional
    title : str, optional
    figsize : tuple of float, float
    adjust : bool, optional
    compress_dups : bool, optional
    """
    if not folder_name:
        folder_name = data.name

    utils.make_folder(folder_name)

    if not title:
        title = "{} - {} - (Bio-N={}, Tech-N={})".format(
            data.tissue,
            data.enrichment,
            min(len(group) for group in data.groups.values()),
            data.sets,
        )
        if abs(fold_cutoff - 1.2) > 0.1:
            title += " -- Fold Change > {}".format(fold_cutoff)

    if title:
        file_name = re.sub("[ ></]", "_", title) + "_Volcano.png"

        if folder_name:
            file_name = os.path.join(folder_name, file_name)

    upper_fold = np.log2(fold_cutoff)
    lower_fold = -upper_fold

    # Calculate the Fold-Change / p-values
    pvals = []
    changes = []
    sig_pvals = []
    sig_changes = []
    sig_labels = []
    colors = []

    for index, row in data.psms.iterrows():
        color = "grey"
        row_pval = -np.log10(row["p-value"])
        row_change = np.log2(row["Fold Change"])

        pvals.append(row_pval)
        changes.append(row_change)

        if row_pval > pval_cutoff and \
           (row_change > upper_fold or row_change < lower_fold):
            row_label = " / ".join(sorted(row["Proteins"].genes))

            sig_pvals.append(row_pval)
            sig_changes.append(row_change)
            sig_labels.append(row_label)
            color = "blue"

        colors.append(color)

    if compress_dups:
        sig_pvals, sig_changes, sig_labels = _remove_lesser_dups(
            sig_pvals, sig_changes, sig_labels
        )

    # Draw the figure
    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(changes, pvals, c=colors)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("$log_2$ Fold Change")
    ax.set_ylabel("$-log_{10}$ p-value")
    ax.axhline(pval_cutoff, color="r", linestyle="dashed", linewidth=0.5)
    ax.axvline(upper_fold, color="r", linestyle="dashed", linewidth=0.5)
    ax.axvline(lower_fold, color="r", linestyle="dashed", linewidth=0.5)
    ax.set_ylim(bottom=-0.1)

    # Position the labels
    texts = []
    for x, y, txt in zip(sig_changes, sig_pvals, sig_labels):
        text = ax.text(x, y, txt)

        if highlight and txt in highlight:
            text.set_fontsize(20)

        text.set_bbox(
            dict(
                color="lightgreen" if x > 0 else "pink",
                alpha=0.8,
                edgecolor="red",
            )
        )

        texts.append(text)

    if adjust:
        adjust_text(
            x=sig_changes,
            y=sig_pvals,
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

    if title:
        ax.set_title(title)
        fig.savefig(
            file_name,
            bbox_inches="tight", dpi=500,
            transparent=True,
        )

    fig.show()


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

    f.show()

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

    f.show()

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
                    snr_cutoff=1,
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
    data, sequences,
):
    """
    Plot the levels of a sequence between two groups.

    Parameters
    ----------
    data : :class:`DataSet<pyproteome.data_sets.DataSet>`
    sequences : list of str
    """
    groups = list(reversed(data.groups.values()))
    labels = list(reversed(data.groups.keys()))

    psms = data.psms.copy()
    psms["Seq Str"] = psms["Sequence"].apply(str)
    psms = psms[psms["Seq Str"].isin(sequences)]

    points = np.array(
        [
            psms[group].as_matrix().sum(axis=0)
            for group in groups
        ]
    )

    values = points.mean(axis=1)
    errs = points.std(axis=1)

    f, ax = plt.subplots()

    indices = np.arange(len(values))
    bar_width = .35
    ax.bar(
        bar_width + indices,
        values,
        bar_width,
        yerr=errs,
        ecolor="k",
    )

    ax.set_ylabel(
        "Cumulative Channel Signal{}".format(
            " (Normalized)" if data.normalized else ""
        ),
        fontsize=20,
    )
    ax.ticklabel_format(style="plain")

    for label in ax.get_yticklabels():
        label.set_fontsize(14)

    ax.set_xticks(indices + bar_width * 1.5)
    ax.set_xticklabels(labels, fontsize=16)

    title = "{}".format(
        " / ".join(sequences),
    )
    ax.set_title(title, fontsize=20)
    ax.xaxis.grid(False)

    def _wrap_list(val):
        if isinstance(val, float):
            return [val]
        return val

    # XXX: Not correct...
    display(
        dict(
            zip(
                sequences,
                _wrap_list(
                    ttest_ind(
                        points[0, :].T,
                        points[1, :].T
                    )[1]
                )
            )
        )
    )

    return f


def plot_sequence(
    data, sequence,
):
    """
    Plot the levels of a sequence across multiple channels.

    Parameters
    ----------
    data : :class:`DataSet<pyproteome.data_sets.DataSet>`
    sequence : str or :class:`Sequence<pyproteome.sequence.Sequence>`
    """
    channels = list(data.channels.keys())

    psms = data.psms[data.psms["Sequence"] == sequence]

    values = psms[channels].as_matrix()
    values = (values.T / values[:, 0]).T

    f, ax = plt.subplots()

    for i in range(values.shape[0]):
        indices = np.arange(len(values[i]))
        bar_width = .35
        ax.bar(bar_width + indices, values[i], bar_width)
        ax.set_xticks(indices + bar_width * 1.5)
        ax.set_xticklabels(list(data.channels.values()))

    display(values)


def plot_correlation(
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

        text = ax.text(
            x, y, " / ".join(row["Proteins_x"].genes),
        )

        text.set_bbox(
            dict(
                color="lightgreen" if ratio < 0 else "pink",
                alpha=0.8,
                edgecolor="red",
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

    ax.set_xlabel("$log_2$ Fold Change -- {}".format(data1.name))
    ax.set_ylabel("$log_2$ Fold Change -- {}".format(data2.name))

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

    tfs = data.psms[data.psms["Proteins"].apply(_is_tf)]
    tfs.sort(columns="Fold Change", ascending=False, inplace=True)

    if csv_name:
        tfs[["Proteins", "Sequence", "Modifications", "Fold Change"]].to_csv(
            csv_name,
            index=False,
        )

    return tfs
