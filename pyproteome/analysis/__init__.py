# -*- coding: UTF-8 -*-
"""
This module provides functionality for data set analysis.

Functions include volcano plots, sorted tables, and plotting sequence levels.
"""

from __future__ import division

# Built-ins
import logging
import os

# IPython
# from IPython.display import display

# Core data analysis libraries
# from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns


from . import fetch_data, utils


LOGGER = logging.getLogger("pyproteome.analysis")


def hierarchical_heatmap(
    data,
    baseline_channels=None,
    metric="euclidean", method="centroid",
    row_cluster=True, col_cluster=True
):
    """
    Plot a hierarhically-clustered heatmap of a data set.

    Parameters
    ----------
    data : :class:`DataSet<pyproteome.data_sets.DataSet>`
    baseline_channels : list of str, optional
        List of channels to average and use as baseline for each row.
    metric : str, optional
        Hierarchical clustering distance metric.
    method : str, optional
        Hierarchical clustering method.
    row_cluster : bool, optional
    col_cluster : bool, optional
    """
    channels = [
        channel
        for channels in data.groups.values()
        for channel in channels
    ]
    channel_names = [
        data.channels[i]
        for i in channels
    ]

    if baseline_channels is None:
        baseline_channels = list(data.groups.values())[-1]

    psms = data.psms[
        data.psms["Fold Change"].apply(lambda x: max([x, 1/x])) > 1.5
    ]

    raw = psms[channels].as_matrix()
    raw_baseline = psms[baseline_channels].as_matrix().mean(axis=1)
    raw = np.log2((raw.T / raw_baseline).T)

    sns.clustermap(
        raw,
        method=method,
        metric=metric,
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        xticklabels=channel_names,
        yticklabels=False,
    )


def find_tfs(data, folder_name=None, csv_name=None):
    """
    Scan over a data set to find proteins annotated as transcription factors.

    Parameters
    ----------
    data : :class:`DataSet<pyproteome.data_sets.DataSet>`
    folder_name : str, optional
    csv_name : str, optional
    """
    if folder_name is None:
        folder_name = data.name

    utils.make_folder(folder_name)

    if csv_name is None:
        csv_name = "Changing TFs.csv"

    if folder_name and csv_name:
        csv_name = os.path.join(folder_name, csv_name)

    def _is_tf(prots):
        go_terms = (
            "DNA binding",
            "double-stranded DNA binding",
            "transcription factor binding",
            "transcription, DNA-templated",
        )
        return any(
            go_term in go
            for prot in prots
            for go in fetch_data.get_uniprot_data(prot.accession).get(
                "go", [],
            )
            for go_term in go_terms
        )

    tfs = data.psms[data.psms["Proteins"].apply(_is_tf)].copy()
    tfs.sort_values(by="Fold Change", ascending=False, inplace=True)

    if csv_name:
        tfs[["Proteins", "Sequence", "Modifications", "Fold Change"]].to_csv(
            csv_name,
            index=False,
        )

    tfs = tfs[["Proteins", "Sequence", "Fold Change", "p-value"]]

    return tfs.style.set_table_styles(  # Hide index and "Validated" columns
        [
            {"selector": "th:first-child", "props": [("display", "none")]},
            {"selector": "*", "props": [("text-align", "left")]},
        ]
    )
