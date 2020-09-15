
from __future__ import division

import pickle
import logging

import numpy as np
import pandas as pd

import brainrnaseq as brs
from . import cache

LOGGER = logging.getLogger('brainrnaseq.enrichments')


def build_barres_table(cell_types=None, force=False):
    if not force:
        try:
            with open(cache.ENRICHMENT_CACHE, 'rb') as f:
                return pickle.load(f)
        except:
            pass

    if cell_types is None:
        cell_types = brs.DEFAULT_CELL_TYPES

    cache.get_barres_seq_data(force=force)

    LOGGER.info('Calculating cell type enrichment')
    enriched = {}

    for species, data in cache.BARRES_SPECIES_DATA.items():
        enriched[species] = {}

        for col in data.columns:
            if hasattr(col, 'lower') and col.lower() == 'gene':
                gene_col_name = col
                break

        data = data.drop_duplicates(subset=gene_col_name)

        for _, row in data.iterrows():
            means = {
                cell_type: np.nanmean(np.array(row[cols].values, dtype=float))
                for cell_type, cols in brs.CELL_TYPE_COLS[species].items()
                if cell_type in cell_types
            }
            max_cell, max_mean = max(means.items(), key=lambda x: x[1])
            enrichment = max_mean / (sum(means.values()) - max_mean)

            enriched[species][row[gene_col_name]] = max_cell, enrichment

    try:
        with open(cache.ENRICHMENT_CACHE, 'wb') as f:
            pickle.dump(enriched, f)
    except Exception as e:
        LOGGER.warning(
            'Unable to save enrichment information to cache: {}'
            .format(e)
        )

    return enriched


def _fix_name(name):
    name = name.title()
    return {
        'Endothelial': 'Endothelia',
        'Endothelial Cells': 'Endothelia',
        'Myeloid': 'Microglia',
        'Oligodendrocyte': 'Myelinating Oligodendrocytes',
        'Oligodendrocyte Precursor Cells': 'OPC',
    }.get(name, name)


def build_hansen_table(
    cell_types=None,
    force=False,
    column=None,
):
    if cell_types is None:
        cell_types = brs.DEFAULT_CELL_TYPES

    cache.get_hansen_seq_data(force=force)

    LOGGER.info('Calculating cell type enrichment')
    enriched = {}

    for species, data in cache.HANSEN_SPECIES_DATA.items():
        if column is None:
            col = {
                'Homo sapiens': 'Barres Human Cell Types',
                'Mus musculus': 'Barres Mouse Cell Types',
            }.get(species)
        else:
            col = column

        enriched[species] = {
            ind: (_fix_name(row[col]), 10)
            for ind, row in data.iterrows()
            if not pd.isnull(row[col])
        }

    return enriched


def get_enrichments(
    species,
    add_mappings=True,
    cutoff=2.5,
    backend='Barres',
    **kwargs
):
    '''
    Fetch enrichment table mapping genes to cell types with enriched expression
    in the brain.

    Parameters
    ----------
    species : str
    add_mappings : bool, optional
    cutoff : float, optional
    backend : str, optional
        One of {'Barres', 'Hansen'}.

    Returns
    -------
    dict of (str, str)
    '''
    build_table = {
        'Barres': build_barres_table,
        'Hansen': build_hansen_table,
    }.get(backend)

    enrichments = {
        key: val[0]
        for key, val in build_table(**kwargs)[species].items()
        if val[1] >= cutoff
    }

    if add_mappings:
        for index, row in cache.get_mapping_data(species=species).iterrows():
            if index not in enrichments:
                continue

            for syn in row['Synonyms'].split('|'):
                if syn in enrichments or syn == '-':
                    continue

                enrichments[syn] = enrichments[index]

    return enrichments
