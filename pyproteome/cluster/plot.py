
from __future__ import absolute_import, division

from functools import cmp_to_key
import logging
from math import ceil, sqrt

from fastcluster import linkage
from matplotlib import pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram
import sklearn

import pyproteome as pyp

from . import clusterer

LOGGER = logging.getLogger('pyp.cluster.plot')
COLOR_MAP = plt.cm.rainbow
CORR_COLOR_MAP = plt.cm.coolwarm


def hierarchical_clusters(data, y_pred):
    ds = data['ds']
    z = data['z']

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
            z[y_pred.values == x].mean(axis=0)
            for x in ss
        ])
    )

    leaves = dendrogram(
        link,
        labels=[str(x) for x in ss],
        no_plot=True,
    )['leaves']

    mapping = {
        val: key
        for key, val in enumerate(leaves)
    }

    new_ind = np.argsort(np.vectorize(lambda x: mapping[x])(y_pred))
    cn = np.corrcoef(data['data'].values[new_ind])

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
    div_scale=None,
    show_names=None,
):
    if ax is None:
        f, ax = plt.subplots(figsize=(13, 12))
    elif f is None:
        f = ax.get_figure()

    if div_scale is None:
        div_scale = max([data['data'].shape[0] // 1000, 1])

    if len(set(y_pred)) > 1:
        mapping, new_ind, _, cn = hierarchical_clusters(data, y_pred)
        sorted_clusters = sorted(mapping, key=lambda x: mapping[x])
        y_pred = y_pred.values[new_ind]
    else:
        sorted_clusters = sorted(set(y_pred))
        cn = data['c']

    if div_scale > 1:
        cn = cn[::div_scale, ::div_scale]

    mesh = ax.pcolormesh(
        cn,
        cmap=CORR_COLOR_MAP,
        vmin=-1, vmax=1,
    )

    if show_names:
        sorted_clusters = [
            i
            for i in sorted_clusters
            if i in show_names
        ]

    for cluster_n in sorted_clusters:
        ind = np.arange(
            0,
            data['data'].shape[0],
        )[y_pred == cluster_n]
        xy = np.median(ind) / div_scale

        if show_names and isinstance(show_names, dict):
            cluster_n = show_names.get(cluster_n)

        ax.text(
            s=str(cluster_n),
            x=xy,
            y=xy,
            # fontsize=sqrt(ind.shape[0]),
            color='k',
            ha='center',
            va='center',
        )

    if colorbar:
        f.colorbar(mesh)

    return f


def plot_cluster(
    data, y_pred, cluster_n,
    f=None,
    ax=None,
    div_scale=None,
    ylabel=True,
    title=None,
    color=None,
):
    if ax is None:
        f, ax = plt.subplots(
            dpi=pyp.DEFAULT_DPI,
            figsize=(8, 4),
        )
    elif f is None:
        f = ax.get_figure()

    z = data['z']
    names = data['names']

    dad = z[y_pred.values == cluster_n]
    n = dad.shape[0]

    if div_scale is None:
        div_scale = max([dad.shape[0] // 100, 1])

    if div_scale > 1:
        dad = dad[::div_scale]

    prev_ind = 0
    means = np.mean(dad, axis=0)

    for ind, v in enumerate(data['classes']):
        n_v, p_v = None, None

        if ind > 0:
            p_v = data['classes'][ind - 1]

        if ind + 1 < len(data['classes']):
            n_v = data['classes'][ind + 1]

        if v != p_v and p_v is not None:
            ax.axvline(
                x=ind - .5,
                color='k',
                linestyle=':',
            )

        if v == n_v:
            continue

        ax.plot(
            np.tile(
                np.arange(prev_ind, ind + 1),
                (dad.shape[0], 1),
            ).T,
            dad[:, prev_ind:ind + 1].T,
            alpha=min([5 / dad.shape[0], 1]),
            color=color or COLOR_MAP(cluster_n / len(set(y_pred))),
        )

        ax.plot(
            range(prev_ind, ind + 1),
            means[prev_ind:ind + 1],
            color='k',
        )

        prev_ind = ind + 1

    if ylabel:
        ax.set_ylabel('Z-scored Change')

    ax.set_facecolor((.9,) * 3)
    ax.set_xticks(range(dad.shape[1]))
    ax.set_xticklabels(names, rotation=45, horizontalalignment='right')
    ax.set_title(
        'Cluster #{}, N={}'.format(cluster_n, n)
        if title is None else
        title
    )

    return ax


def plot_all_clusters(
    data, y_pred,
    cols=4,
):
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
        )

    for ax in ax_iter:
        ax.axis('off')

    f.tight_layout(h_pad=1)

    return f


def show_cluster(
    data, y_pred,
    seq=None,
    protein=None,
    ylabel=True,
    f=None,
    ax=None,
    color=None,
    div_scale=None,
):
    if ax is None:
        f, ax = plt.subplots(figsize=(6, 6))
    elif f is None:
        f = ax.get_figure()

    ds = data['ds']
    z = data['z']
    dp = ds.copy()
    mod = ''
    mask = None

    if seq is not None:
        mask = dp['Sequence'] == seq
        mod = dp[mask].iloc[0]['Modifications']
        protein = ' / '.join(dp.filter(series=mask).genes)
        cluster = y_pred[mask].iloc[0]
    elif protein is not None:
        mask = dp['Proteins'] == protein
        mod = dp[mask].iloc[0]['Modifications']
        protein = ' / '.join(dp.filter(series=mask).genes)
        cluster = y_pred[mask].iloc[0]

    mod_str = mod.get_mods([(None, 'Phospho')]).__str__(prot_index=0)

    title = '{}{}'.format(
        (protein[:15] + ' ...') if len(protein) > 15 else protein,
        ' {}'.format(mod_str) if mod_str else '',
    )
    plot_cluster(
        data, y_pred, cluster,
        ax=ax,
        color=color,
        div_scale=div_scale,
        ylabel=ylabel,
        title=title,
    )

    if mask is not None:
        z_m = z[mask]
        prev_ind = 0

        for ind, v in enumerate(data['classes']):
            n_v = None

            if ind + 1 < len(data['classes']):
                n_v = data['classes'][ind + 1]

            if v == n_v:
                continue

            ax.plot(
                np.tile(
                    np.arange(prev_ind, ind + 1),
                    (z_m.shape[0], 1),
                ).T,
                z_m[:, prev_ind:ind + 1].T,
                color='k',
                linestyle=':',
                linewidth=5,
            )

            prev_ind = ind + 1

    dp.psms = dp[y_pred == cluster]

    return dp, f, ax


def show_peptide_clusters(
    data, y_pred,
    filters,
    new_colors=False,
    div_scale=None,
    cols=4,
):
    rows = int(ceil(len(filters) / cols))
    f, axes = plt.subplots(
        rows, cols,
        figsize=(cols * 2, rows * 2),
        sharex=True,
        sharey=True,
    )
    ax_iter = iter(axes.ravel())

    clusters = [
        sorted(
            set(
                y_pred[
                    data['ds']['Sequence'] == fil['seq']
                    if 'seq' in fil else
                    data['ds']['Proteins'] == fil['protein']
                ]
            )
        )
        for fil in filters
    ]

    clusters = [
        i[0]
        for i in clusters
        if i
    ]

    ss = sorted(set(clusters))

    colors = [
        COLOR_MAP(ss.index(cluster_n) / len(ss))
        if new_colors else
        COLOR_MAP(cluster_n / len(set(y_pred)))
        for cluster_n in clusters
    ]

    for ind, (fil, color, ax) in enumerate(zip(filters, colors, ax_iter)):
        show_cluster(
            data, y_pred,
            ax=ax,
            ylabel=ind % cols == 0,
            color=color,
            div_scale=div_scale,
            save=False,
            **fil
        )

    for ax in ax_iter:
        ax.axis('off')

    return f, axes


def pca(data):
    classes = data['classes']
    z = data['z']
    names = data['names']

    f, ax = plt.subplots(1, 1, figsize=(4, 4))
    x = sklearn.decomposition.PCA().fit_transform(z.T)

    for ind, label in enumerate(data['labels']):
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
    ax.set_title('PCA {}'.format(data['ds'].name))
    ax.set_xlabel('Component 1')
    ax.set_ylabel('Component 2')

    return f, ax


def cluster_range(
    data,
    min_clusters=2, max_clusters=20,
    cols=3,
):
    clusters = range(min_clusters, max_clusters + 1)

    rows = int(np.ceil(len(clusters) / cols))

    f, axes = plt.subplots(
        rows,
        cols,
        figsize=(4 * cols, 4 * rows),
    )

    for n, ax in zip(clusters, axes.ravel()):
        LOGGER.info('Generating plots for cluster {}'.format(n))

        _, y_pred = clusterer.cluster(
            data,
            n_clusters=n,
        )

        pyp.cluster.plot.cluster_corrmap(
            data, y_pred,
            ax=ax,
            colorbar=False,
        )

        ax.set_title('{} Clusters'.format(n))

    for ax in axes.ravel()[len(clusters):]:
        ax.axis('off')

    return f, axes
