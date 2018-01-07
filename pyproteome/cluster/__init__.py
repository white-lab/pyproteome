
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

    f, ax = plt.subplots(1, 1, figsize=(6, 6))
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
