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
PSP_SITE_MAPPING_URL = (
    "https://www.phosphosite.org/downloads/Phosphorylation_site_dataset.gz"
)
ORGANISM_MAPPING = {
    # 'cat': ,
    # 'chicken': ,
    "Bos taurus": 'cow',
    "Canis familiaris": 'dog',
    "Mustela putorius": 'ferret',
    # 'frog': ,
    "Drosophila melanogaster": 'fruit fly',
    # 'goat': ,
    # 'guinea pig': ,
    # 'hamster': ,
    "Equus caballus": 'horse',
    "Homo sapiens": 'human',
    # 'monkey': ,
    "Mus musculus": 'mouse',
    # 'papillomavirus': ,
    # 'pig': ,
    # 'quail': ,
    # 'rabbit': ,
    "Rattus norvegicus": 'rat',
    # 'sheep': ,
    # 'starfish': ,
    # 'torpedo': ,
    # 'turkey': ,
    # 'water buffalo': ,
}


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


@pyp.utils.memoize
def _get_phosphomap_data():
    url = PSP_SITE_MAPPING_URL

    r = requests.get(url, stream=True)
    r.raise_for_status()

    gz = gzip.GzipFile(fileobj=io.BytesIO(r.raw.read()))

    return pd.read_table(gz, skiprows=[0, 1, 2], sep="\t")


def _remap_data(ds, species):
    new = ds.copy()
    old_species = list(new.species)[0]

    old_species = ORGANISM_MAPPING.get(old_species, old_species)
    species = ORGANISM_MAPPING.get(species, species)

    new.species = set([species])

    LOGGER.info(
        "Remapping phosphosites from {} to {}".format(old_species, species)
    )

    mapping = _get_phosphomap_data()

    acc_mapping = mapping.set_index(["ACC_ID", "MOD_RSD", "ORGANISM"])
    site_mapping = mapping.set_index(["SITE_GRP_ID", "ORGANISM"])

    acc_mapping = acc_mapping.sort_index()
    site_mapping = site_mapping.sort_index()

    def _remap(val):
        acc, mod = val.split(",")

        try:
            site = acc_mapping.loc[acc, mod, old_species]
        except KeyError:
            return None

        site = site.iloc[0]["SITE_GRP_ID"]

        try:
            re_map = site_mapping.loc[site, species]
        except KeyError:
            return None

        re_map = re_map.iloc[0]

        return ",".join([re_map["ACC_ID"], re_map["MOD_RSD"]])

    new.psms["ID"] = new.psms["ID"].apply(_remap)
    new.psms = new.psms[~(new.psms["ID"].isnull())]

    LOGGER.info(
        "Mapping {} IDs ({}) down to {} IDs ({})"
        .format(ds.shape[0], old_species, new.shape[0], species)
    )

    return new


def _get_psite_ids(ds, species):
    LOGGER.info("Building list of individual phosphosites")
    new_rows = []

    def _split_rows(row):
        acc = row["Proteins"].accessions[0]

        for mod in row["Modifications"].get_mods("Phospho"):
            new_row = row.to_dict()

            mod_name = "{}{}-p".format(mod.letter, mod.abs_pos[0] + 1)

            new_row["ID"] = ",".join([acc, mod_name])
            new_rows.append(new_row)

    ds.psms.apply(_split_rows, axis=1)

    return pd.DataFrame(new_rows)


def _get_protein_ids(ds, species):
    return ds.psms["Proteins"].apply(
        lambda row:
        brs.mapping.get_entrez_mapping(
            row.genes[0],
            species=species,
        ),
    )


def _get_pathways(species):
    LOGGER.info("Fetching WikiPathways")

    pathways_df = _get_wikipathways(species)

    if species in ["Homo sapiens"]:
        pathways_df = pathways_df.append(_get_pathway_common(species))

    if species in ["Mus musculus"]:
        pathways_df = pathways_df.append(_get_gskb_pathways(species))

    return pathways_df


def gsea(
    ds, phenotype,
    p_sites=False,
    remap=None, folder_name=None,
    **kwargs
):
    folder_name = pyp.utils.make_folder(
        data=ds,
        folder_name=folder_name,
        sub="Volcano",
    )

    LOGGER.info(
        "filtering ambiguous peptides {}".format(len(set(ds.genes)))
    )

    ds = ds.filter(fn=lambda x: len(x["Proteins"].genes) == 1)
    species = list(ds.species)[0]

    LOGGER.info("building gene sets")
    pathways_df = _get_pathways(species)
    LOGGER.info("Loaded {} gene sets".format(pathways_df.shape[0]))

    LOGGER.info("mapping genes: {}, {}".format(len(set(ds.genes)), p_sites))

    if p_sites:
        ds.psms = _get_psite_ids(ds, species)

        if remap:
            ds = _remap_data(ds, remap)
    else:
        ds.psms["ID"] = _get_protein_ids(ds, species)

    LOGGER.info("building correlations")

    phenotype = pd.to_numeric(phenotype)

    ds.psms = ds.psms[
        (~ds.psms["ID"].isnull()) &
        (
            (~ds.psms[phenotype.index].isnull()).sum(axis=1) >=
            enrichments.MIN_PERIODS
        )
    ]

    ds = enrichments.correlate_phenotype(ds, phenotype)

    LOGGER.info("plotting enrichments {}".format(ds.shape))

    vals = enrichments.plot_enrichment(
        ds, pathways_df, phenotype,
        **kwargs
    )

    vals.to_csv(
        os.path.join(
            folder_name,
            ("PSEA" if p_sites else "GSEA") + ".csv",
        ),
    )

    return vals
