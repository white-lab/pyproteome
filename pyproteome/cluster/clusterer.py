
from scipy.stats import zscore

import numpy as np
import pandas as pd
import sklearn
import sklearn.decomposition
import sklearn.cluster

import pyproteome as pyp


def get_data(ds, dropna=True, corrcoef=True, groups=None):
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

    if corrcoef:
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
        "c": c if corrcoef else None,
        "names": names,
        "labels": groups,
        "classes": classes,
    }


def cluster(data, z=False, log2=False, clr=None, n_clusters=20):
    if z:
        x = data["z"]
    else:
        if log2:
            x = data['ds'].data.applymap(np.log2)
        else:
            x = data['ds'].data

    if clr is None:
        clr = sklearn.cluster.AgglomerativeClustering(n_clusters=n_clusters)

    y_pred_spec = clr.fit_predict(x)

    y_pred = pd.Series(
        [
            i if i >= 0 else max(y_pred_spec) + 1
            for i in y_pred_spec
        ],
        index=data["ds"].psms.index,
    )

    return clr, y_pred
