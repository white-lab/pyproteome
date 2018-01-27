# -*- coding: UTF-8 -*-
"""
This module provides functionality for data set analysis.

Functions include volcano plots, sorted tables, and plotting sequence levels.
"""

from __future__ import absolute_import, division

# Built-ins
import logging
import os

# IPython
# from IPython.display import display

# Core data analysis libraries
# from matplotlib import pyplot as plt


import pyproteome as pyp


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
