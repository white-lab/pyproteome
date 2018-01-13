
from __future__ import division

from math import ceil, sqrt
from functools import cmp_to_key
import os

from pyproteome import utils

from matplotlib import pyplot as plt
import numpy as np

from scipy.cluster.hierarchy import dendrogram
from fastcluster import linkage


def hierarchical_clusters(data, y_pred):
    ds = data["ds"]
    z = data["z"]

    ss = sorted(
        set(y_pred),
        key=cmp_to_key
        (
            lambda x, y:
            np.corrcoef(
                ds[y_pred == x].mean(axis=0),
                ds[y_pred == y].mean(axis=0),
            )[0, 1]
        )
    )

    mapping = {
        val: key
        for key, val in enumerate(ss)
    }

    link = linkage(
        np.array([
            z[y_pred.as_matrix() == x].mean(axis=0)
            for x in sorted(set(y_pred))
        ])
    )

    leaves = dendrogram(
        link,
        labels=[str(x) for x in sorted(set(y_pred))],
        no_plot=True,
    )["leaves"]

    mapping = {
        val: key
        for key, val in enumerate(leaves)
    }

    new_ind = np.argsort(np.vectorize(lambda x: mapping[x])(y_pred))
    cn = np.corrcoef(data["data"].as_matrix()[new_ind])

    cmapping = {
        cluster_n: index
        for index, cluster_n in enumerate(
            sorted(mapping, key=lambda x: mapping[x])
        )
    }

    return mapping, new_ind, cmapping, cn


def _get_folder(data):
    folder_name = os.path.join(
        data["ds"].name, "Clusters",
    )
    utils.make_folder(folder_name)
    return folder_name


def cluster_corrmap(
    data, y_pred,
    colorbar=True,
    f=None,
    ax=None,
):
    folder_name = _get_folder(data)
    if ax is None:
        f, ax = plt.subplots(figsize=(13, 12))
    elif f is None:
        f = ax.get_figure()

    div_scale = max([data["data"].shape[0] // 1000, 5])

    mapping, new_ind, cmapping, cn = hierarchical_clusters(data, y_pred)

    mesh = ax.pcolormesh(
        cn[::div_scale, ::div_scale],
        cmap=plt.cm.Spectral_r,
        vmin=-1, vmax=1,
    )

    for cluster_n in sorted(mapping, key=lambda x: mapping[x]):
        ind = np.arange(
            0,
            data["data"].shape[0],
        )[y_pred.as_matrix()[new_ind] == cluster_n]
        xy = np.median(ind) / div_scale
        ax.text(
            s=str(cluster_n),
            x=xy,
            y=xy,
            fontsize=sqrt(ind.shape[0]),
            color='k',
            horizontalalignment='center',
            verticalalignment='center',
        )

    if colorbar:
        f.colorbar(mesh)

    f.savefig(
        os.path.join(folder_name, "Cluster-Corrmap.png"),
        bbox_inches="tight", dpi=300,
        transparent=True,
    )
    return f


def plot_cluster(
    data, y_pred, cluster_n,
    f=None,
    ax=None,
    ylabel=True,
    title=None,
    color=None,
):
    if ax is None:
        f, ax = plt.subplots()
    elif f is None:
        f = ax.get_figure()

    z = data["z"]
    names = data["names"]

    dad = z[y_pred.as_matrix() == cluster_n]
    n = dad.shape[0]
    div_scale = max([dad.shape[0] // 100, 1])
    dad = dad[::div_scale]

    for i in range(dad.shape[0]):
        ax.plot(
            dad[i, :],
            alpha=min([1 / dad.shape[0] * 50 / 2, 1]),
            color=color or plt.cm.rainbow(cluster_n / max(y_pred)),
        )

    ax.plot(
        np.mean(dad, axis=0),
        color='k',
    )

    for ind, (v, n_v) in enumerate(zip(data["classes"], data["classes"][1:])):
        if v == n_v:
            continue

        ax.axvline(
            x=ind + .5,
            color="k",
            linestyle="--",
        )

    if ylabel:
        ax.set_ylabel("Z-scored Change")

    ax.set_xticks(range(dad.shape[1]))
    ax.set_xticklabels(names, rotation=45, horizontalalignment="right")
    ax.set_title(
        "Cluster #{}, N={}".format(cluster_n, n)
        if title is None else
        title
    )


def plot_all_clusters(
    data, y_pred,
    cols=4,
):
    folder_name = _get_folder(data)
    clusters = sorted(set(y_pred))
    f, axes = plt.subplots(
        int(np.ceil(len(clusters) / cols)),
        cols,
        figsize=(cols * 3, int(np.ceil(len(clusters) / cols)) * 3),
        sharey=True,
        sharex=True,
    )

    ax_iter = iter(axes.ravel())

    for ind, (cluster_n, ax) in enumerate(zip(sorted(set(y_pred)), ax_iter)):
        plot_cluster(
            data, y_pred, cluster_n,
            ax=ax,
            ylabel=ind % cols == 0,
        )

    for ax in ax_iter:
        ax.axis('off')

    f.tight_layout(h_pad=1)

    f.savefig(
        os.path.join(folder_name, "Clusters.png"),
        bbox_inches="tight", dpi=300,
        transparent=True,
    )


def show_cluster(
    data, y_pred,
    seq=None,
    protein=None,
    ylabel=True,
    f=None,
    ax=None,
    color=None,
):
    folder_name = _get_folder(data)

    if ax is None:
        f, ax = plt.subplots(figsize=(6, 6))
    elif f is None:
        f = ax.get_figure()

    ds = data["ds"]
    z = data["z"]
    dp = ds.copy()
    mod = ""

    if seq is not None:
        cluster = list(set([i for i in y_pred[dp["Sequence"] == seq]]))[0]
        mod = dp[dp["Sequence"] == seq].iloc[0]["Modifications"]
        protein = " / ".join(
            gene
            for prot in dp[dp["Sequence"] == seq]["Proteins"]
            for gene in prot.genes
        )
    elif protein is not None:
        cluster = list(set([i for i in y_pred[dp["Proteins"] == protein]]))[0]

    title = "Cluster #{} N={}\n({}{})".format(
        cluster,
        (dp[y_pred == cluster]).shape[0],
        protein,
        " {}".format(mod.__str__(prot_index=0)) if mod else "",
    )
    plot_cluster(
        data, y_pred, cluster,
        ax=ax,
        color=color,
        ylabel=ylabel,
        title=title,
    )

    if seq is not None:
        ax.plot(
            z[(dp["Sequence"] == seq).as_matrix()][0],
            color="k",
            linestyle="--",
            linewidth=5,
        )
    elif protein is not None:
        ax.plot(
            z[(dp["Proteins"] == protein).as_matrix()][0],
            color="k",
            linestyle="--",
            linewidth=5,
        )

    f.savefig(
        os.path.join(folder_name, "Cluster-{}.png".format(cluster)),
        bbox_inches="tight", dpi=300,
        transparent=True,
    )

    dp.psms = dp[y_pred == cluster]

    return dp


def show_peptide_clusters(
    data, y_pred,
    filters,
    new_colors=True,
    cols=4,
):
    folder_name = _get_folder(data)

    rows = int(ceil(len(filters) / cols))
    f, axes = plt.subplots(
        rows, cols,
        figsize=(cols * 4, rows * 4),
        sharex=True,
        sharey=True,
    )
    ax_iter = iter(axes.ravel())

    clusters = [
        sorted(
            set(
                y_pred[
                    data["ds"]["Sequence"] == fil["seq"]
                    if "seq" in fil else
                    data["ds"]["Proteins"] == fil["protein"]
                ]
            )
        )[0]
        for fil in filters
    ]

    ss = sorted(set(clusters))

    colors = [
        plt.cm.rainbow(ss.index(cluster_n) / len(ss))
        if new_colors else
        plt.cm.rainbow(cluster_n / max(y_pred))
        for cluster_n in clusters
    ]

    for ind, (fil, color, ax) in enumerate(zip(filters, colors, ax_iter)):
        show_cluster(
            data, y_pred,
            ax=ax,
            ylabel=ind % cols == 0,
            color=color,
            **fil
        )

    for ax in ax_iter:
        ax.axis('off')

    f.savefig(
        os.path.join(folder_name, "PeptideClusters.png"),
        bbox_inches="tight", dpi=300,
        transparent=True,
    )

    return f, axes
