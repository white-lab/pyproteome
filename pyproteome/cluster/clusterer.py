
from scipy.stats import zscore

import numpy as np
import pandas as pd
import sklearn
import sklearn.decomposition
import sklearn.cluster

import pyproteome as pyp


def get_data(ds, dropna=True, groups=None):
    ds = ds.copy()

    if groups is None:
        groups = ds.cmp_groups or [list(ds.groups.keys())]

    if dropna:
        ds = ds.dropna(how="any", groups=pyp.utils.flatten_list(groups))

    names = [
        chan
        for lst in groups
        for group in lst
        for chan in ds.groups[group]
        if chan in ds.channels
    ]
    chans = [ds.channels[chan] for chan in names]
    data = ds.data[chans]

    c = np.corrcoef(data.values)
    z = zscore(data, axis=1)

    classes = np.array([
        [
            lst.index(i)
            for lst in groups
            for i in lst
            if col in ds.groups[i]
        ][0]
        for col in names
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


def cluster(data, clr=None, n_clusters=20):
    z = data["z"]
    ds = data["ds"]

    if clr is None:
        clr = sklearn.cluster.AgglomerativeClustering(n_clusters=n_clusters)

    y_pred_spec = clr.fit_predict(z)

    y_pred = pd.Series(
        [
            i if i >= 0 else max(y_pred_spec) + 1
            for i in y_pred_spec
        ],
        index=ds.psms.index,
    )

    return clr, y_pred


def cluster_clusters(
    data, y_pred,
    p_cutoff=1,
    corr_cutoff=.7,
    max_iterations=50,
):
    y_pred = y_pred.copy()

    for _ in range(max_iterations):
        ss = sorted(set(y_pred), key=lambda i: -(y_pred == i).sum())
        renamed = False

        for ind, i in enumerate(ss):
            corrs = [
                (
                    o,
                    [
                        np.correlate(
                            data["z"][y_pred == i].mean(axis=0),
                            data["z"][y_pred == o].mean(axis=0),
                        ),
                        0,
                    ],
                )
                for o in sorted(set(y_pred))
                if o < i
            ]

            if not corrs:
                continue

            o, (rho, p) = max(corrs, key=lambda x: x[1][0])

            if (
                p < p_cutoff and
                rho > 0 and rho > corr_cutoff
            ):
                y_pred[y_pred == i] = o
                renamed = True

        if not renamed:
            break

    y_pred_old = y_pred.copy()

    for ind, i in enumerate(sorted(set(y_pred))):
        y_pred[y_pred_old == i] = ind

    return y_pred
