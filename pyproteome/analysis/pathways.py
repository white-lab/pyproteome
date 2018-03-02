# -*- coding: UTF-8 -*-
"""
This module provides functionality for data set analysis.

Functions include volcano plots, sorted tables, and plotting sequence levels.
"""

from __future__ import absolute_import, division

# Built-ins
import logging
import os
import requests

# IPython
# from IPython.display import display

# Core data analysis libraries
# from matplotlib import pyplot as plt
import numpy as np
import pandas as pd


import pyproteome as pyp
import brainrnaseq as brs
from . import correlation, enrichments


LOGGER = logging.getLogger("pyproteome.pathways")


def find_tfs(data, folder_name=None, csv_name=None):
    """
    Scan over a data set to find proteins annotated as transcription factors.

    Parameters
    ----------
    data : :class:`DataSet<pyproteome.data_sets.DataSet>`
    folder_name : str, optional
    csv_name : str, optional
    """
    folder_name = pyp.utils.make_folder(
        data=data,
        folder_name=folder_name,
        sub="Pathway Analysis",
    )

    if csv_name is None:
        csv_name = "Changing TFs.csv"

    uni_data = pyp.fetch_data.fetch_uniprot_data(data.accessions)

    def _is_tf(row):
        go_terms = (
            "DNA binding",
            "double-stranded DNA binding",
            "transcription factor binding",
            "transcription, DNA-templated",
        )
        return any(
            go_term in go
            for acc in row["Proteins"].accessions
            for go in uni_data.get(acc, {}).get("go", [])
            for go_term in go_terms
        )

    tfs = data.filter(fn=_is_tf).psms
    tfs.sort_values(by="Fold Change", ascending=False, inplace=True)

    if csv_name:
        tfs[["Proteins", "Sequence", "Modifications", "Fold Change"]].to_csv(
            os.path.join(folder_name, csv_name),
            index=False,
        )

    tfs = tfs[["Proteins", "Sequence", "Fold Change", "p-value"]]

    return tfs.style.set_table_styles(  # Hide index and "Validated" columns
        [
            {"selector": "th:first-child", "props": [("display", "none")]},
            {"selector": "*", "props": [("text-align", "left")]},
        ]
    )


def _get_pathways():
    url = (
        "https://raw.githubusercontent.com/dhimmel/pathways/"
        "master/data/pathways.tsv"
    )
    # response = requests.get(url, stream=True)
    # response.raise_for_status()
    #
    # with open("pathways.tsv", "wb") as f:
    #     for block in response.iter_content(1024):
    #         f.write(block)
    #
    pathways_df = pd.read_table(url)

    return pathways_df


def gsea(ds, phenotype, **kwargs):
    LOGGER.info("filtering ambiguous peptides {}".format(len(set(ds.genes))))
    ds = ds.copy()
    ds = ds.filter(fn=lambda x: len(x["Proteins"].genes) == 1)

    LOGGER.info("building gene sets")
    pathways_df = _get_pathways()

    pathways_df["set"] = pathways_df["genes"].apply(
        lambda row:
        set(int(i) for i in row.split("|")),
    )

    LOGGER.info("mapping genes: {}".format(len(set(ds.genes))))

    ds.psms["Entrez"] = ds.psms["Proteins"].apply(
        lambda row:
        brs.mapping.get_entrez_mapping(
            row.genes[0],
            species="Human",
        ),
    )

    LOGGER.info("building correlations")

    phenotype = pd.to_numeric(phenotype)

    ds.psms["Correlation"] = ds.psms.apply(
        lambda row:
        phenotype.corr(
            pd.to_numeric(row[phenotype.index]),
            method="pearson",
            min_periods=5,
        ),
        # correlation.spearmanr_nan(
        #     row[phenotype.index],
        #     phenotype,
        # ).correlation,
        axis=1,
    )

    ds.psms = ds.psms[
        (~ds.psms["Entrez"].isnull()) &
        ((~ds.psms[phenotype.index].isnull()).sum(axis=1) >= 5)
    ]

    LOGGER.info("plotting enrichments {}".format(ds.shape))

    return enrichments.plot_enrichment(
        ds, pathways_df, phenotype,
        **kwargs
    )
