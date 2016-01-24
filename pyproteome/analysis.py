"""
This module provides functionality for data set analysis.

Functions include volcano plots, sorted tables, and plotting sequence levels.
"""

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
import networkx as nx
import numpy as np
# import pandas as pd
# import seaborn as sns
# import scipy
from scipy.stats import ttest_ind
# from scipy.stats.mstats import mquantiles
# from scipy.cluster import hierarchy
# import sklearn
# from sklearn.cluster import KMeans, MiniBatchKMeans

# Misc extras
from adjustText import adjust_text
import Bio.Alphabet.IUPAC
import Bio.Seq
import Bio.motifs
# import fastcluster as fst
# import somoclu
# import uniprot

from . import utils, fetch_data, modifications, sequence


LOGGER = logging.getLogger("pyproteome.analysis")


def snr_table(
    psms,
    snr_cutoff=None, fold_cutoff=None,
    folder_name=None, csv_name=None,
):
    """
    Show a signal to noise table.

    Parameters
    ----------
    psms : pandas.DataFrame
    snr_cutoff : float, optional
    fold_cutoff : float, optional
    folder_name : str, optional
    csv_name : str, optional
    """
    utils.make_folder(folder_name)

    if folder_name and csv_name:
        csv_name = os.path.join(folder_name, csv_name)

    psms = psms[["Proteins", "Sequence", "SNR", "Fold Change"]]

    # psms["Sort"] = psms["SNR"].apply(abs)
    psms["Sort"] = psms["Fold Change"].apply(lambda x: max([x, 1/x]))

    psms.sort_values("Sort", inplace=True, ascending=False)
    psms.drop("Sort", axis=1, inplace=True)

    if csv_name:
        psms.to_csv(csv_name)

    if snr_cutoff:
        psms = psms[abs(psms["SNR"]) > snr_cutoff]

    if fold_cutoff:
        psms = psms[
            np.maximum.reduce(
                [
                    psms["Fold Change"],
                    1 / psms["Fold Change"],
                ]
            ) > fold_cutoff
        ]

    display(psms)


def _place_labels(x, y, texts, ax=None, spring_k=None, spring_scale=None):
    if spring_k is None:
        spring_k = 0.15

    if spring_scale is None:
        spring_scale = 0.1

    if ax is None:
        ax = plt.gca()

    if not texts:
        return

    G = nx.Graph()
    init_pos = {}
    data_pos = {}
    data_nodes = []
    ano_nodes = []

    for j, b in enumerate(zip(x, y, texts)):
        x, y, label = b
        data_str = "data_{}".format(j)
        ano_str = "ano_{}".format(j)

        G.add_node(ano_str)
        G.add_node(data_str)
        G.add_edge(ano_str, data_str, weight=0.1)

        data_nodes.append(data_str)
        ano_nodes.append(ano_str)

        data_pos[data_str] = (x, y)
        init_pos[data_str] = (x, y)
        init_pos[ano_str] = (x, y)

    pos = nx.spring_layout(
        G, pos=init_pos, fixed=data_nodes,
        k=spring_k,
        scale=spring_scale,
    )

    for j, txt in enumerate(texts):
        data_str = "data_{}".format(j)
        ano_str = "ano_{}".format(j)

        ax.annotate(
            txt,
            xy=data_pos[data_str], xycoords="data",
            xytext=pos[ano_str], textcoords="data",
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            # fontsize=20,
            bbox=dict(boxstyle='square', fc='pink', ec='none'),
        )


def volcano_plot(
    norm_psms,
    pval_cutoff=1.3, fold_cutoff=1.2, folder_name=None, title=None,
    spring_k=None, spring_scale=None, adjust_layout=True,
):
    """
    Display a volcano plot of data.

    This plot inclues the fold-changes and p-values associated with said
    changes.

    Parameters
    ----------
    norm_psms : pandas.DataFrame
    pval_cutoff : float, optional
    fold_cutoff : float, optional
    folder_name : str, optional
    title : str, optional
    spring_k : float, optional
    spring_scale : float, optional
    adjust_layout : bool, optional
        Use the adjustText library to position labels, otherwise use networkx
        and its spring_layout function.
    """
    utils.make_folder(folder_name)

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

    for index, row in norm_psms.iterrows():
        color = "grey"
        row_pval = -np.log10(row["p-value"])
        row_change = np.log2(row["Fold Change"])

        pvals.append(row_pval)
        changes.append(row_change)

        if row_pval > pval_cutoff and \
           (row_change > upper_fold or row_change < lower_fold):
            row_label = " / ".join(row["Proteins"].genes)

            sig_pvals.append(row_pval)
            sig_changes.append(row_change)
            sig_labels.append(row_label)
            color = "blue"

        colors.append(color)

    # Draw the figure
    fig, ax = plt.subplots()  # figsize=(12, 12))
    ax.scatter(changes, pvals, c=colors)
    ax.set_xlabel("$log_2$ Fold Change")
    ax.set_ylabel("$-log_{10}$ p-value")
    ax.axhline(pval_cutoff, color="r", linestyle="dashed", linewidth=0.5)
    ax.axvline(upper_fold, color="r", linestyle="dashed", linewidth=0.5)
    ax.axvline(lower_fold, color="r", linestyle="dashed", linewidth=0.5)
    ax.set_ylim(bottom=-0.1)

    # Position the labels
    if adjust_layout:
        adjust_text(
            x=sig_changes,
            y=sig_pvals,
            texts=sig_labels,
            ax=ax,
        )
    else:
        _place_labels(
            x=sig_changes,
            y=sig_pvals,
            texts=sig_labels,
            ax=ax,
            spring_k=spring_k,
            spring_scale=spring_scale,

        )

    if title:
        ax.set_title(title)
        fig.savefig(
            file_name,
            bbox_inches='tight', dpi=500,
            transparent=True,
        )

    fig.show()


def venn3(psms_a, psms_b, psms_c, folder_name=None, filename=None):
    """
    Display a three-way venn diagram.

    Parameters
    ----------
    psms_a : pandas.DataFrame
    psms_b : pandas.DataFrame
    psms_c : pandas.DataFrame
    folder_name : str, optional
    filename : str, optional
    """
    utils.make_folder(folder_name)

    if folder_name and filename:
        filename = os.path.join(folder_name, filename)

    group_a = set(psms_a)
    group_b = set(psms_b)
    group_c = set(psms_c)

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
        set_labels=("Hippocampus", "Cortex", "Cerebellum"),
    )

    for label in v.set_labels:
        if label:
            label.set_fontsize(32)

    for label in v.subset_labels:
        if label:
            label.set_fontsize(20)

    f.show()

    if filename:
        f.savefig(filename, transparent=True)


def significant_sequences(psms, snr_cutoff=1, fold_cutoff=None, p_cutoff=None):
    """
    Process psms and keep those with significant changes.

    Removes any row with a SNR, Fold Change, or p-value above or below a given
    threshold.

    Parameters
    ----------
    psms : pandas.DataFrame
    snr_cutoff : float, optional
    fold_cutoff : float, optional
    p_cutoff : float, optional

    Returns
    -------
    pandas.DataFrame
    """
    if fold_cutoff is None:
        fold_cutoff = 1.3

        if snr_cutoff < 0:
            fold_cutoff = 1 / fold_cutoff

    if snr_cutoff < 0:
        psms = psms[psms["SNR"] < snr_cutoff]
    else:
        psms = psms[psms["SNR"] > snr_cutoff]

    if p_cutoff:
        psms = psms[psms["p-value"] < p_cutoff]

    if fold_cutoff < 1:
        psms = psms[psms["Fold Change"] < fold_cutoff]
    else:
        psms = psms[psms["Fold Change"] > fold_cutoff]

    return psms["Sequence"]


def write_lists(
    psms,
    folder_name=None, sorted_name="sorted_list.txt",
    hits_name="hits_list.txt", background_name="back_list.txt",
):
    """
    Write a list of peptides to files.

    Includes peptides sorted by fold change, significantly changing peptides,
    and background peptides with little change.

    Parameters
    ----------
    psms : pandas.DataFrame
    folder_name : str, optional
    sorted_name : str, optional
    hits_name : str, optional
    background_name : str, optional
    """
    utils.make_folder(folder_name)

    if folder_name:
        sorted_name = os.path.join(folder_name, sorted_name)
        hits_name = os.path.join(folder_name, hits_name)
        background_name = os.path.join(folder_name, background_name)

    change_psms = psms.copy()
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
                for i in psms[
                    np.logical_and(
                        abs(psms["SNR"]) >= 1,
                        np.maximum.reduce(
                            [
                                change_psms["Fold Change"],
                                1 / change_psms["Fold Change"],
                            ]
                        ) > 1.3
                    )
                ]["Proteins"].drop_duplicates(keep="first")
            )
        )

    with open(background_name, "w") as f:
        f.write(
            "\n".join(
                i.accessions[0]
                for i in psms[
                    abs(psms["SNR"]) < 0.1
                ]["Proteins"].drop_duplicates(keep="first")
            )
        )


def plot_sequence_between(
    psms, sequences,
    group_a, group_b,
    xlabels=None,
    normalize=True,
):
    """
    Plot the levels of a sequence between two groups.

    Parameters
    ----------
    psms : pandas.DataFrame
    sequences : list of str
    group_a : list of str
    group_b : list of str
    normalize : bool, optional
    """
    if normalize:
        group_a = utils.norm(group_a)
        group_b = utils.norm(group_b)

    psms = psms.copy()
    psms["Seq Str"] = psms["Sequence"].apply(str)
    psms = psms[psms["Seq Str"].isin(sequences)]

    values_a = psms[group_a].as_matrix()
    values_b = psms[group_b].as_matrix()

    values = [values_a.mean(), values_b.mean()]
    errs = [values_a.std(), values_b.std()]

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

    ax.set_ylabel("Cumulative TMT Signal", fontsize=20)
    ax.ticklabel_format(style="plain")

    for label in ax.get_yticklabels():
        label.set_fontsize(14)

    if xlabels:
        ax.set_xticks(indices + bar_width * 1.5)
        ax.set_xticklabels(xlabels, fontsize=16)

    title = "{}".format(
        " / ".join(sequences),
    )
    ax.set_title(title, fontsize=20)
    ax.xaxis.grid(False)

    def _wrap_list(val):
        if isinstance(val, float):
            return [val]
        return val

    display(
        dict(zip(sequences, _wrap_list(ttest_ind(values_a.T, values_b.T)[1])))
    )

    return f


def plot_sequence(
    psms, sequence, channels,
    normalize=True,
):
    """
    Plot the levels of a sequence across multiple channels.

    Parameters
    ----------
    psms : pandas.DataFrame
    sequence : str or pyproteome.Sequence
    channels : list of str
    normalize : bool, optional
    """
    if normalize:
        channels = utils.norm(channels)

    psms = psms[psms["Sequence"] == sequence]

    values = psms[list(channels.keys())].as_matrix()
    values = (values.T / values[:, 0]).T

    f, ax = plt.subplots()

    for i in range(values.shape[0]):
        indices = np.arange(len(values[i]))
        bar_width = .35
        f.bar(bar_width + indices, values[i], bar_width)
        ax.set_xticks(indices + bar_width * 1.5)
        ax.set_xticklabels(list(channels.values()))

    display(values)


def find_tfs(psms, folder_name=None, csv_name=None):
    """
    Scan over psms to find proteins annotated as transcription factors.

    Parameters
    ----------
    psms : pandas.DataFrame
    folder_name : str, optional
    csv_name : str, optional
    """
    if folder_name and csv_name:
        csv_name = os.path.join(folder_name, csv_name)

    def _is_tf(prots):
        go_terms = (
            "DNA binding",
            "double-stranded DNA binding",
            "transcription factor binding",
        )
        return any(
            go_term in go
            for prot in prots
            for go in fetch_data.get_uniprot_data(prot.accession).get(
                "go", [],
            )
            for go_term in go_terms
        )

    tfs = psms[psms["Proteins"].apply(_is_tf)]
    tfs.sort(columns="Fold Change", ascending=False, inplace=True)

    if csv_name:
        tfs[["Proteins", "Sequence", "Modifications", "Fold Change"]].to_csv(
            csv_name,
            index=False,
        )

    return tfs


def motif_analysis(
    psms, letter_mod_types=None,
    folder_name=None, filename="Motif.svg",
):
    """
    Create a sequence logo figure.

    Logos are created based on the frequencies of peptides in a data set.

    Parameters
    ----------
    psms : pandas.DataFrame
    letter_mod_types : list of tuple of str, str
    folder_name : str, optional
    filename : str, optional
    """
    if folder_name and filename:
        filename = os.path.join(folder_name, filename)

    alpha = Bio.Alphabet.IUPAC.protein
    m = Bio.motifs.create(
        [
            Bio.Seq.Seq(seq.upper(), alphabet=alpha)
            for seq in sequence.generate_n_mers(
                psms["Sequence"],
                letter_mod_types=letter_mod_types,
            )
        ],
        alphabet=alpha,
    )
    m.weblogo(filename, format="svg", logo_title="Testing")
