# -*- coding: UTF-8 -*-
"""
This module provides functionality for data set analysis.

Functions include volcano plots, sorted tables, and plotting sequence levels.
"""

from __future__ import division

# Built-ins
from collections import Iterable, OrderedDict
import itertools
import logging
import os
import re

# IPython
# from IPython.display import display

# Core data analysis libraries
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr, spearmanr, ttest_ind

# Misc extras
from adjustText.adjustText import adjust_text

from . import fetch_data, utils


LOGGER = logging.getLogger("pyproteome.analysis")
DEFAULT_DPI = 300


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
    figsize=(12, 8),
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
    psms["Sort Ind"] = psms["Seq Str"].apply(lambda x: sequences.index(x))
    psms = psms.sort_values("Sort Ind")

    values = np.array([
        (
            np.nanmean(psms[channel].as_matrix(), axis=1) /
            np.nanmean(psms[channels[0]].as_matrix(), axis=1)
        )
        for channel in channels
    ])
    means = values

    errs = [
        np.nanstd(psms[channel].as_matrix(), axis=1) /
        np.nanmean(psms[channels[0]].as_matrix(), axis=1)
        for channel in channels
    ]
    if psms.shape[0] <= 1:
        means = [[i] for i in values]
        errs = [[i] for i in errs]

    values, labels = zip(*[
        (value, label)
        for label, value in zip(data.groups.keys(), values)
        if value.any() and
        (cmp_groups is None or any(label in i for i in cmp_groups))
    ])

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

    f, ax = plt.subplots(figsize=figsize)

    indices = np.arange(len(means))
    width = .75
    bar_width = width / len(means[0])

    for ind in range(len(means[0])):
        ax.bar(
            bar_width * ind - width / 2 + indices,
            [i[ind] for i in means],
            bar_width,
            # yerr=[i[ind] for i in errs],
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

    ax.set_xticks(indices - bar_width / 2)
    ax.set_xticklabels(labels, fontsize=20)

    title = "{}".format(
        # " / ".join(sequences),
        " / ".join(
            sorted(
                set(
                    gene
                    for row in psms["Proteins"]
                    for gene in row.genes
                )
            )
        )
    )[:50]
    ax.set_title(title, fontsize=32)
    ax.xaxis.grid(False)

    # if not means[0].any():
    #     return f

    if len(means[0]) > 1:
        pep_seqs = [
            data["Sequence"][data["Sequence"] == seq].iloc[0]
            for seq in sequences
        ]
        ax.legend([
            "{} ({} - {})".format(
                str(seq),
                seq.protein_matches[0].rel_pos,
                seq.protein_matches[0].rel_pos + len(seq.pep_seq),
            )
            for seq in pep_seqs
        ])
        return f, ax

    y_max = np.max([max(i) + max(j) for i, j in zip(means, errs)])

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

    return f, ax


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
    channel_names = list(OrderedDict.fromkeys(channel_names))
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
    ).dropna(subset=("Fold Change_x", "Fold Change_y"))

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

    name1 = data1.name
    name2 = data2.name

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
        f.savefig(
            filename,
            transparent=True,
            dpi=DEFAULT_DPI,
        )


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
    p_cutoff=0.05,
    fold_cutoff=1.25,
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

    if title is None:
        title = data.name

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

        if row["corr p-value"] < p_cutoff:
            sig_x.append(row["Correlation"])
            sig_y.append(-np.log10(row["corr p-value"]))
            sig_labels.append(" / ".join(sorted(row["Proteins"].genes)))

        colors.append(
            "blue"
            if row["corr p-value"] < p_cutoff else
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

    if title:
        ax.set_title(title, fontsize=32)

    cp.psms = cp.psms[cp.psms["corr p-value"] < p_cutoff]

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

    if not isinstance(datas, Iterable):
        datas = [datas]

    if protein:
        seqs = list(
            set(
                str(seq)
                for data in datas
                for seq in data[data["Proteins"] == protein]["Sequence"]
            )
        )
        datas = [
            i.filter(proteins=[protein])
            for i in datas
        ]

    if isinstance(seqs, str):
        seqs = [seqs]

    for seq in seqs:
        for data in datas:
            try:
                os.makedirs(data.name)
            except:
                pass

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
                    ),
                    figsize=figsize,
                )
                f.savefig(
                    os.path.join(
                        data.name,
                        re.sub(
                            "[?/]",
                            "_",
                            "{}-{}.png".format(
                                prot,
                                data.name,
                            ),
                        ),
                    ),
                    bbox_inches="tight",
                    dpi=DEFAULT_DPI,
                    transparent=True,
                )

            if between:
                f, _ = plot_sequence_between(
                    data, [seq],
                    cmp_groups=cmp_groups,
                )
                f.savefig(
                    os.path.join(
                        data.name,
                        re.sub(
                            "[?/]",
                            "_",
                            "{}-{}-between.png".format(
                                prot,
                                data.name,
                            ),
                        ),
                    ),
                    bbox_inches="tight",
                    dpi=DEFAULT_DPI,
                    transparent=True,
                )


def plot_all_together(
    datas, seqs=None, protein=None, proteins=None, only=True, **kwargs
):
    assert seqs is not None or protein is not None or proteins is not None

    if protein:
        seqs = list(
            set(
                str(seq)
                for data in datas
                for seq in data[
                    data["Proteins"].apply(
                        lambda x:
                        all(gene == protein for gene in x.genes)
                    )
                    if only else
                    data["Proteins"] == protein
                ]["Sequence"]
            )
        )
    elif proteins:
        seqs = list(
            set(
                str(seq)
                for data in datas
                for seq in data[data["Proteins"].apply(
                    lambda x: any(x == i for i in proteins)
                )]["Sequence"]
            )
        )
        print(seqs)

    if isinstance(seqs, str):
        seqs = [seqs]

    for data in datas:
        try:
            os.makedirs(data.name)
        except:
            pass

        prot = " / ".join(
            gene
            for i in data[data["Sequence"] == seqs[0]]["Proteins"]
            for gene in i.genes
        )
        f, ax = plot_sequence_between(
            data, seqs,
            **kwargs
        )
        f.savefig(
            os.path.join(
                data.name,
                re.sub(
                    "[?/]",
                    "_",
                    "{}-{}-between.png".format(
                        prot,
                        data.name,
                    ),
                ),
            ),
            bbox_inches="tight",
            dpi=DEFAULT_DPI,
            transparent=True,
        )

    return f, ax
