
from __future__ import absolute_import, division

from math import ceil, sqrt
from functools import cmp_to_key
import os

from fastcluster import linkage
from matplotlib import pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram
from scipy.stats import zscore
import seaborn as sns
import sklearn

import pyproteome as pyp

from . import clusterer


def _make_folder(data=None, folder_name=None):
    if folder_name is None:
        folder_name = os.path.join(
            data.name
            if data is not None else
            os.getcwd(),
            "Clusters",
        )

    return pyp.utils.makedirs(folder_name)


def hierarchical_heatmap(
    data,
    cmp_groups=None,
    baseline_channels=None,
    metric="euclidean",
    method="centroid",
    row_cluster=True,
    col_cluster=True,
    show_y=False,
    filename="Hierarchical Heatmap.png",
    folder_name=None,
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
    folder_name = _make_folder(data, folder_name=folder_name)

    if cmp_groups is None:
        cmp_groups = [list(data.groups.keys())]

    channels = [
        data.channels[channel]
        for groups in cmp_groups
        for group in groups
        for channel in data.groups[group]
        if channel in data.channels
    ]

    raw = data.psms[channels].T.apply(zscore).T
    raw.index = data.psms["Proteins"].apply(str)
    raw = raw.dropna(how="all")

    map = sns.clustermap(
        raw,
        method=method,
        metric=metric,
        row_cluster=row_cluster,
        col_cluster=col_cluster,
    )
    map.ax_heatmap.set_xticklabels(
        map.ax_heatmap.get_xticklabels(),
        rotation=45,
    )

    if not show_y:
        map.ax_heatmap.set_yticklabels([])
        map.ax_heatmap.set_ylabel("")

    if filename:
        map.savefig(
            os.path.join(folder_name, filename),
            bbox_inches="tight",
            dpi=pyp.DEFAULT_DPI,
            transparent=True,
        )

    return map


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


def cluster_corrmap(
    data, y_pred,
    colorbar=True,
    f=None,
    ax=None,
    filename="Cluster-Corrmap.png",
    folder_name=None,
):
    folder_name = _make_folder(data["ds"], folder_name=folder_name)

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

    if filename:
        f.savefig(
            os.path.join(folder_name, filename),
            bbox_inches="tight",
            dpi=pyp.DEFAULT_DPI,
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
    folder_name=None,
):
    folder_name = _make_folder(data["ds"], folder_name=folder_name)

    ss = sorted(set(y_pred))
    rows = int(np.ceil(len(ss) / cols))
    f, axes = plt.subplots(
        rows, cols,
        figsize=(cols * 3, rows * 3),
        sharey=True,
        sharex=True,
    )

    ax_iter = iter(axes.ravel())

    for ind, (cluster_n, ax) in enumerate(zip(ss, ax_iter)):
        plot_cluster(
            data, y_pred, cluster_n,
            ax=ax,
            ylabel=ind % cols == 0,
            save=False,
        )

    for ax in ax_iter:
        ax.axis('off')

    f.tight_layout(h_pad=1)

    f.savefig(
        os.path.join(folder_name, "Clusters.png"),
        bbox_inches="tight",
        dpi=pyp.DEFAULT_DPI,
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
    save=True,
    folder_name=None,
):
    folder_name = _make_folder(data["ds"], folder_name=folder_name)

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

    if save:
        f.savefig(
            os.path.join(folder_name, "Cluster-{}.png".format(cluster)),
            bbox_inches="tight",
            dpi=pyp.DEFAULT_DPI,
            transparent=True,
        )

    dp.psms = dp[y_pred == cluster]

    return dp


def show_peptide_clusters(
    data, y_pred,
    filters,
    new_colors=True,
    cols=4,
    folder_name=None,
):
    folder_name = _make_folder(data["ds"], folder_name=folder_name)

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
            save=False,
            **fil
        )

    for ax in ax_iter:
        ax.axis('off')

    f.savefig(
        os.path.join(folder_name, "PeptideClusters.png"),
        bbox_inches="tight",
        dpi=pyp.DEFAULT_DPI,
        transparent=True,
    )

    return f, axes


def pca(data):
    classes = data["classes"]
    z = data["z"]
    names = data["names"]
    print(z)

    f, ax = plt.subplots(1, 1, figsize=(4, 4))
    x = sklearn.decomposition.PCA().fit_transform(z.T)

    for ind, label in enumerate(data["labels"]):
        ax.scatter(
            x[classes == ind, 0],
            x[classes == ind, 1],
            label=label,
        )

    offset = 0
    for ind, column in enumerate(names):
        ax.text(
            x=x[ind, 0] + offset,
            y=x[ind, 1] + offset,
            s=column,
        )

    ax.legend()
    ax.set_title("PCA {}".format(data["ds"].name))
    ax.set_xlabel("Component 1")
    ax.set_ylabel("Component 2")


def cluster_range(
    data,
    min_clusters=2, max_clusters=20,
    cols=3,
    folder_name=None,
    filename="Cluster-Range-Scan.png",
):
    folder_name = _make_folder(data["ds"], folder_name=folder_name)

    clusters = range(min_clusters, max_clusters + 1)

    rows = int(np.ceil(len(clusters) / cols))

    f, axes = plt.subplots(
        rows,
        cols,
        figsize=(4 * cols, 4 * rows),
    )

    for n, ax in zip(clusters, axes.ravel()):
        print(n)

        _, y_pred = clusterer.cluster(
            data,
            n_clusters=n,
        )

        pyp.cluster.plot.cluster_corrmap(
            data, y_pred,
            ax=ax,
            colorbar=False,
        )

        ax.set_title("{} Clusters".format(n))

    for ax in axes.ravel()[len(clusters):]:
        ax.axis("off")

    f.savefig(
        os.path.join(folder_name, filename),
        bbox_inches="tight",
        dpi=pyp.DEFAULT_DPI,
        transparent=True,
    )
