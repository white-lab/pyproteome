
import os

import pyproteome
from pyproteome import cluster, volcano, motifs, tables


def auto_clusterer(
    data,
    get_data_kwargs=None,
    cluster_kwargs=None,
    cluster_cluster_kwargs=None,
    folder_name="Clusters",
):
    """
    Cluster and generate plots for a data set.

    Parameters
    ----------
    data : :class:`DataSet<pyproteome.data_sets.DataSet>`
    """
    get_data_kwargs = get_data_kwargs or {}
    cluster_kwargs = cluster_kwargs or {}
    cluster_cluster_kwargs = cluster_cluster_kwargs or {}

    data = cluster.get_data(
        data,
        **get_data_kwargs
    )

    cluster.pca(data)

    if "n_clusters" not in cluster_kwargs:
        cluster_kwargs["n_clusters"] = 100

    _, y_pred_old = cluster.cluster(
        data,
        **cluster_kwargs
    )

    y_pred = cluster.cluster_clusters(
        data, y_pred_old,
        **cluster_cluster_kwargs
    )

    cluster.plot.cluster_corrmap(data, y_pred_old)
    cluster.plot.cluster_corrmap(data, y_pred)

    cluster.plot.plot_all_clusters(data, y_pred)

    ss = sorted(set(y_pred))

    for ind in ss:
        volcano.plot_volcano_filtered(
            data["ds"], {"series": y_pred == ind},
            title="Cluster {}".format(ind),
            folder_name=folder_name,
        )
        f = motifs.logo.make_logo(
            data["ds"], {"series": y_pred == ind},
            title="Cluster {}".format(ind),
        )
        f.savefig(
            os.path.join(folder_name, "Logo - Cluster {}.png".format(ind)),
            bbox_inches="tight",
            dpi=pyproteome.DEFAULT_DPI,
            transparent=True,
        )

    slices = [
        data["ds"].filter({"series": y_pred == ind})
        for ind in ss
    ]

    for ind, s in zip(ss, slices):
        s.name = "Cluster {}".format(ind)

    tables.write_full_tables(
        slices,
        folder_name="Clusters",
    )
