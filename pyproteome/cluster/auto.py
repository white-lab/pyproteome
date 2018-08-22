
from __future__ import absolute_import, division

import os
import pickle

from matplotlib import pyplot as plt

import pyproteome as pyp


def auto_clusterer(
    data,
    get_data_kwargs=None,
    cluster_kwargs=None,
    cluster_cluster_kwargs=None,
    plot_clusters_kwargs=None,
    volcano_kwargs=None,
    plots=True,
    close=False,
    folder_name=None,
    filename="clusters.pkl",
):
    """
    Cluster and generate plots for a data set.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.DataSet`
    """
    folder_name = pyp.utils.make_folder(
        data=data,
        folder_name=folder_name,
        sub="Clusters",
    )

    get_data_kwargs = get_data_kwargs or {}
    cluster_kwargs = cluster_kwargs or {}
    cluster_cluster_kwargs = cluster_cluster_kwargs or {}
    plot_clusters_kwargs = plot_clusters_kwargs or {}
    volcano_kwargs = volcano_kwargs or {}

    data = pyp.cluster.get_data(
        data,
        **get_data_kwargs
    )

    if "n_clusters" not in cluster_kwargs:
        cluster_kwargs["n_clusters"] = 100

    clr, y_pred_old = pyp.cluster.cluster(
        data,
        **cluster_kwargs
    )

    y_pred = pyp.cluster.cluster_clusters(
        data, y_pred_old,
        **cluster_cluster_kwargs
    )

    if not plots:
        return data, y_pred

    pyp.cluster.plot.pca(data)

    pyp.cluster.plot.cluster_corrmap(
        data, y_pred_old,
        filename="First Clusters.png",
    )
    pyp.cluster.plot.cluster_corrmap(
        data, y_pred,
        filename="Final Clusters.png",
    )

    pyp.cluster.plot.plot_all_clusters(
        data, y_pred,
        **plot_clusters_kwargs
    )

    # return data, y_pred, clr

    ss = sorted(set(y_pred))

    for ind in ss:
        f = pyp.cluster.plot.plot_cluster(
            data, y_pred, ind,
            div_scale=1,
        )
        if f:
            f.savefig(
                os.path.join(folder_name, "Cluster-{}.png".format(ind)),
                bbox_inches="tight",
                dpi=pyp.DEFAULT_DPI,
                transparent=True,
            )
            if close:
                plt.close(f)

        f, _ = pyp.volcano.plot_volcano_filtered(
            data["ds"], {"series": y_pred == ind},
            title="Cluster {}".format(ind),
            folder_name=folder_name,
            **volcano_kwargs
        )

        if f and close:
            plt.close(f)

        f, _ = pyp.motifs.logo.make_logo(
            data["ds"], {"series": y_pred == ind},
            title="Cluster {}".format(ind),
        )

        if f:
            f.savefig(
                os.path.join(folder_name, "Logo - Cluster {}.png".format(ind)),
                bbox_inches="tight",
                dpi=pyp.DEFAULT_DPI,
                transparent=True,
            )
            if close:
                plt.close(f)

    slices = [
        data["ds"].filter({"series": y_pred == ind})
        for ind in ss
    ]

    for ind, s in zip(ss, slices):
        s.name = "Cluster {}".format(ind)

    pyp.tables.write_full_tables(
        slices,
        folder_name=folder_name,
    )

    if filename:
        with open(os.path.join(folder_name, filename), "wb") as f:
            pickle.dump(y_pred, f)

    return data, y_pred, clr
