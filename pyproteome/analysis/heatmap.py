
import os

import numpy as np
import seaborn as sns

import pyproteome as pyp


def hierarchical_heatmap(
    data,
    cmp_groups=None,
    baseline_channels=None,
    metric="euclidean",
    method="centroid",
    row_cluster=True,
    col_cluster=True,
    show_y=False,
    title=None,
    filename="Hierarchical Heatmap.png",
    folder_name=None,
):
    """
    Plot a hierarhically-clustered heatmap of a data set.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.DataSet`
    baseline_channels : list of str, optional
        List of channels to average and use as baseline for each row.
    metric : str, optional
        Hierarchical clustering distance metric.
    method : str, optional
        Hierarchical clustering method.
    row_cluster : bool, optional
    col_cluster : bool, optional

    Returns
    -------
    map : :class:`seaborn.ClusterGrid`
    """
    folder_name = pyp.utils.make_folder(
        data=data,
        folder_name=folder_name,
        sub="Clusters",
    )

    if cmp_groups is None:
        cmp_groups = [list(data.groups.keys())]

    channels = [
        data.channels[channel]
        for groups in cmp_groups
        for group in groups
        for channel in data.groups[group]
        if channel in data.channels
    ]

    raw = data.psms
    raw = raw.sort_values("Fold Change", ascending=False)
    raw = raw[channels]

    def zscore(x):
        return (x - np.nanmean(x)) / np.nanstd(x)

    # raw = raw.apply(zscore, axis=1)
    raw = raw.apply(np.log2, axis=1)
    raw.index = data.psms["Proteins"].apply(str).apply(lambda x: x[:50])
    raw = raw.dropna(how="all")

    minmax = max([
        abs(raw.min().min()),
        abs(raw.max().max()),
    ])
    # minmax = 1.5

    raw = raw.fillna(value=-minmax)

    map = sns.clustermap(
        raw,
        method=method,
        metric=metric,
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        vmin=-minmax,
        vmax=minmax,
    )
    map.ax_heatmap.set_xticklabels(
        map.ax_heatmap.get_xticklabels(),
        rotation=45,
        horizontalalignment="right",
    )

    if title is not None:
        map.fig.suptitle(title)

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
