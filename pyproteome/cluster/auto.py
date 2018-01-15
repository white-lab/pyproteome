
from __future__ import absolute_import, division

import os

import pyproteome as pyp


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
    try:
        os.makedirs(folder_name)
    except:
        pass

    get_data_kwargs = get_data_kwargs or {}
    cluster_kwargs = cluster_kwargs or {}
    cluster_cluster_kwargs = cluster_cluster_kwargs or {}

    data = pyp.cluster.get_data(
        data,
        **get_data_kwargs
    )

    pyp.cluster.pca(data)

    if "n_clusters" not in cluster_kwargs:
        cluster_kwargs["n_clusters"] = 100

    _, y_pred_old = pyp.cluster.cluster(
        data,
        **cluster_kwargs
    )

    y_pred = pyp.cluster.cluster_clusters(
        data, y_pred_old,
        **cluster_cluster_kwargs
    )

    pyp.cluster.plot.cluster_corrmap(data, y_pred_old)
    pyp.cluster.plot.cluster_corrmap(data, y_pred)

    pyp.cluster.plot.plot_all_clusters(data, y_pred)

    ss = sorted(set(y_pred))

    for ind in ss:
        pyp.volcano.plot_volcano_filtered(
            data["ds"], {"series": y_pred == ind},
            title="Cluster {}".format(ind),
            folder_name=folder_name,
        )

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

    slices = [
        data["ds"].filter({"series": y_pred == ind})
        for ind in ss
    ]

    for ind, s in zip(ss, slices):
        s.name = "Cluster {}".format(ind)

    pyp.tables.write_full_tables(
        slices,
        folder_name="Clusters",
    )

    return data, y_pred
