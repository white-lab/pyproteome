
from scipy.stats import zscore

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import sklearn
import sklearn.decomposition
import sklearn.cluster

from . import plot


def get_data(ds, dropna=True, groups=None):
    ds = ds.copy()

    if groups is None:
        groups = list(ds.groups.keys())

    if dropna:
        ds = ds.dropna(how="any", groups=groups)

    data = ds.data
    data = data.loc[:, ~data.columns.duplicated()]
    data = data.loc[
        :,
        data.columns.isin([col for i in groups for col in ds.groups[i]])
    ]

    names = data.columns
    c = np.corrcoef(data.as_matrix())
    z = zscore(data, axis=1)

    classes = np.array([
        [groups.index(i) for i in groups if col in ds.groups[i]][0]
        for col in data.columns
    ])
    print(data.shape)

    return {
        "ds": ds,
        "data": data,
        "z": z,
        "c": c,
        "names": names,
        "labels": groups,
        "classes": classes,
    }


def pca(data):
    classes = data["classes"]
    z = data["z"]
    names = data["names"]

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


def cluster_range(data, min_clusters=2, max_clusters=20, cols=3):
    clusters = range(min_clusters, max_clusters + 1)

    f, axes = plt.subplots(
        int(np.ceil((max_clusters - 1) / cols)),
        cols,
        figsize=(4 * cols, 4 * int(np.ceil((max_clusters - 1) / 2))),
    )
    # axes = np.array([None] * max_clusters)

    for n, ax in zip(clusters, axes.ravel()):
        print(n)

        _, y_pred = cluster(
            data,
            n_clusters=n,
        )

        plot.cluster_corrmap(
            data, y_pred,
            ax=ax,
            colorbar=False,
        )

        ax.set_title("{} Clusters".format(n))

    for ax in axes.ravel()[len(clusters):]:
        ax.axis("off")

    f.savefig(
        "Cluster-Range-Scan.png",
        bbox_inches="tight", dpi=300,
        transparent=True,
    )


def cluster(data, fn=None, kwargs=None, n_clusters=20):
    z = data["z"]
    ds = data["ds"]

    if fn is None:
        fn = sklearn.cluster.AgglomerativeClustering

    if kwargs is None:
        kwargs = {
            "n_clusters": n_clusters,
        }

    clr = fn(
        **kwargs
    )

    y_pred_spec = clr.fit_predict(z)

    y_pred = pd.Series(y_pred_spec, index=ds.psms.index)

    return clr, y_pred


def cluster_clusters(data, y_pred, corr_cutoff=4):
    y_pred = y_pred.copy()

    ss = sorted(set(y_pred))

    for ind, i in enumerate(ss):
        corrs = [
            (
                o,
                np.correlate(
                    data["z"][y_pred == i].mean(axis=0),
                    data["z"][y_pred == o].mean(axis=0),
                )[0],
            )
            for o in sorted(set(y_pred))
            if o < i
        ]
        if not corrs:
            continue

        o, m = max(corrs, key=lambda x: x[1])

        if m > corr_cutoff:
            y_pred[y_pred == i] = o

    y_pred_old = y_pred.copy()

    for ind, i in enumerate(sorted(set(y_pred))):
        y_pred[y_pred_old == i] = ind

    print(sorted(set(y_pred)))
    return y_pred
