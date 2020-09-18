
from __future__ import absolute_import, division

import logging

from matplotlib import pyplot as plt

import pyproteome as pyp

LOGGER = logging.getLogger('pyp.cluster.auto')


def auto_clusterer(
    data,
    get_data_kwargs=None,
    cluster_kwargs=None,
    plot_clusters_kwargs=None,
    volcano_kwargs=None,
    plots=False,
    close=False,
):
    '''
    Cluster and generate plots for a data set.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.DataSet`
    get_data_kwargs : dict
        Arguments passed to :func:`.clusterer.get_data`.
    cluster_kwargs : dict
        Arguments passed to :func:`.clusterer.cluster`.
    plot_clusters_kwargs : dict
        Arguments passed to :func:`.plot.plot_all_clusters`.
    plots : bool, optional
        Generate plots for each cluster.
    close : bool, optional
        Automatically close all figures.

    Returns
    -------
    data : dict
        Dictionary containing the data set and exact matrix used
        for clustering, as well as accessory objects.
    y_pred : :class:`numpy.array`
        List of cluster IDs for each peptide.
    clr
        scikit-learn's cluster object.

    Examples
    --------
    >>> data, y_pred, clr = cluster.auto.auto_clusterer(
    ...     ds,
    ...     get_data_kwargs={
    ...         'dropna': True,
    ...         'corrcoef': False,
    ...     },
    ...     cluster_kwargs={
    ...         'clr': sklearn.cluster.MiniBatchKMeans(
    ...             n_clusters=n,
    ...             random_state=0,
    ...         ),
    ...     },
    ...     plots=False,
    ... )
    '''
    get_data_kwargs = get_data_kwargs or {}
    cluster_kwargs = cluster_kwargs or {}
    plot_clusters_kwargs = plot_clusters_kwargs or {}
    volcano_kwargs = volcano_kwargs or {}

    LOGGER.info('Fetching data matrix.')

    data = pyp.cluster.clusterer.get_data(
        data,
        **get_data_kwargs
    )

    if 'n_clusters' not in cluster_kwargs:
        if 'clr' in cluster_kwargs and hasattr(cluster_kwargs['clr'], 'n_clusters'):
            cluster_kwargs['n_clusters'] = cluster_kwargs['clr'].n_clusters
        else:
            cluster_kwargs['n_clusters'] = 25

    LOGGER.info(
        'Grouping data into n={} clusters.'.format(
            cluster_kwargs['n_clusters'],
        )
    )

    clr, y_pred = pyp.cluster.clusterer.cluster(
        data,
        **cluster_kwargs
    )

    if not plots:
        return data, y_pred, clr

    LOGGER.info(
        'Plotting cluster information (n={} clusters)'.format(len(set(y_pred)))
    )

    # pyp.cluster.plot.pca(data)

    pyp.cluster.plot.cluster_corrmap(
        data, y_pred,
    )

    # pyp.cluster.plot.plot_all_clusters(
    #     data, y_pred,
    #     **plot_clusters_kwargs
    # )

    # # return data, y_pred, clr

    # ss = sorted(set(y_pred))

    # for ind in ss:
    #     LOGGER.info('Plotting cluster #{}'.format(ind))

    #     ax = pyp.cluster.plot.plot_cluster(
    #         data, y_pred, ind,
    #         div_scale=1,
    #     )
    #     f = ax.get_figure()

    #     if f and close:
    #         plt.close(f)
    # #
    # #     f, _ = pyp.volcano.plot_volcano(
    # #         data['ds'].filter(series=y_pred == ind),
    # #         title='Cluster {}'.format(ind),
    # #         **volcano_kwargs
    # #     )[:2]
    # #
    # #     if f and close:
    # #         plt.close(f)
    # #
    # #     f, _ = pyp.motifs.logo.make_logo(
    # #         data['ds'], {'series': y_pred == ind},
    # #         title='Cluster {}'.format(ind),
    # #     )
    # #
    # #     if f and close:
    # #         plt.close(f)

    # slices = [
    #     data['ds'].filter({'series': y_pred == ind})
    #     for ind in ss
    # ]

    # for ind, s in zip(ss, slices):
    #     s.name = 'Cluster {}'.format(ind)

    # pyp.tables.write_full_tables(
    #     slices,
    # )

    return data, y_pred, clr
