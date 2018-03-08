# -*- coding: UTF-8 -*-
"""
This module provides functionality for data set analysis.

Functions include volcano plots, sorted tables, and plotting sequence levels.
"""

from __future__ import absolute_import, division

# Built-ins
import gzip
import io
import logging
import os
import re
import requests

# IPython
# from IPython.display import display

# Core data analysis libraries
# from matplotlib import pyplot as plt
import pandas as pd


import pyproteome as pyp
import brainrnaseq as brs
from . import enrichments


LOGGER = logging.getLogger("pyproteome.pathways")
WIKIPATHWAYS_URL = (
    "http://data.wikipathways.org/20180210/gmt/"
    "wikipathways-20180210-gmt-{}.gmt"
)
PATHWAYS_COMMON_URL = (
    "http://www.pathwaycommons.org/archives/PC2/v9/"
    "PathwayCommons9.All.hgnc.gmt.gz"
)
GSKB_URL = (
    "http://ge-lab.org/gskb/2-MousePath/mGSKB_Entrez.gmt"
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


@pyp.utils.memoize
def _get_gskb_pathways(species):
    LOGGER.info("Fetching GSKB pathways")

    url = GSKB_URL
    r = requests.get(url, stream=True)
    r.raise_for_status()

    def _get_data(line):
        line = line.decode("windows-1252")
        name, _, genes = line.split("\t", 2)

        genes = set(int(i) for i in genes.split("\t") if i)

        return name, genes

    pathways_df = pd.DataFrame(
        data=[
            _get_data(line)
            for ind, line in enumerate(r.iter_lines())
            if ind > 0
        ],
        columns=["name", "set"],
    )

    return pathways_df


@pyp.utils.memoize
def _get_pathway_common(species):
    LOGGER.info("Fetching Pathways Common")

    url = PATHWAYS_COMMON_URL
    r = requests.get(url, stream=True)
    r.raise_for_status()

    name_re = re.compile(
        "name: (.+); datasource: (.+); organism: (.+); idtype: (.+)"
    )

    def _get_data(line):
        line = line.decode()
        _, name, genes = line.split("\t", 2)
        name = name_re.match(name)

        name = {
            "name": name.group(1),
            "datasource": name.group(2),
            "organism": name.group(3),
            "id_type": name.group(4),
        }

        assert int(name["organism"]) == 9606

        genes = set(
            brs.mapping.get_entrez_mapping(i, species=species)
            for i in genes.split("\t")
        )

        return name["name"], genes

    pathways_df = pd.DataFrame(
        data=[
            _get_data(line)
            for line in gzip.GzipFile(fileobj=io.BytesIO(r.raw.read()))
        ],
        columns=["name", "set"],
    )

    return pathways_df


@pyp.utils.memoize
def _get_wikipathways(species):
    url = WIKIPATHWAYS_URL.format("_".join(species.split(" ")))
    response = requests.get(url, stream=True)
    response.raise_for_status()

    def _get_data(line):
        line = line.decode()
        name, _, genes = line.split("\t", 2)
        name, _, _, species = name.split("%")
        assert species == species
        return name, set(int(i) for i in genes.split("\t"))

    pathways_df = pd.DataFrame(
        data=[
            _get_data(line)
            for line in response.iter_lines()
        ],
        columns=["name", "set"],
    )

    return pathways_df


def _get_pathways(species):
    LOGGER.info("Fetching WikiPathways")

    pathways_df = _get_wikipathways(species)

    if species in ["Homo sapiens"]:
        pathways_df = pathways_df.append(_get_pathway_common(species))

    if species in ["Mus musculus"]:
        pathways_df = pathways_df.append(_get_gskb_pathways(species))

    return pathways_df


def gsea(ds, phenotype, **kwargs):
    LOGGER.info("filtering ambiguous peptides {}".format(len(set(ds.genes))))
    ds = ds.filter(fn=lambda x: len(x["Proteins"].genes) == 1)
    species = list(ds.species)[0]

    LOGGER.info("building gene sets")
    pathways_df = _get_pathways(species)
    LOGGER.info("Loaded {} gene sets".format(pathways_df.shape[0]))

    LOGGER.info("mapping genes: {}".format(len(set(ds.genes))))

    ds.psms["Entrez"] = ds.psms["Proteins"].apply(
        lambda row:
        brs.mapping.get_entrez_mapping(
            row.genes[0],
            species=species,
        ),
    )

    LOGGER.info("building correlations")

    phenotype = pd.to_numeric(phenotype)

    ds.psms = ds.psms[
        (~ds.psms["Entrez"].isnull()) &
        (
            (~ds.psms[phenotype.index].isnull()).sum(axis=1) >=
            enrichments.MIN_PERIODS
        )
    ]

    ds = enrichments.correlate_phenotype(ds, phenotype)

    LOGGER.info("plotting enrichments {}".format(ds.shape))

    return enrichments.plot_enrichment(
        ds, pathways_df, phenotype,
        **kwargs
    )
