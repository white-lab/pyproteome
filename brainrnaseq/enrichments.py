
import pickle
import logging

import numpy as np

import brainrnaseq as brs
from . import cache

LOGGER = logging.getLogger("brainrnaseq.enrichments")


def build_enrichment_table(cell_types=None, force=False):
    if not force:
        try:
            with open(cache.ENRICHMENT_CACHE, "rb") as f:
                return pickle.load(f)
        except:
            pass

    if cell_types is None:
        cell_types = brs.DEFAULT_CELL_TYPES

    cache.get_barres_seq_data(force=force)

    LOGGER.info("Calculating cell type enrichment")
    enriched = {}

    for species, data in cache.SPECIES_DATA.items():
        enriched[species] = {}

        for col in data.columns:
            if hasattr(col, "lower") and col.lower() == "gene":
                gene_col_name = col
                break

        data = data.drop_duplicates(subset=gene_col_name)

        for _, row in data.iterrows():
            if any(
                isinstance(row[col], str)
                for cols in brs.CELL_TYPE_COLS[species].values()
                for col in cols
            ):
                continue

            means = {
                cell_type: np.nanmean(np.array(row[cols].values, dtype=float))
                for cell_type, cols in brs.CELL_TYPE_COLS[species].items()
                if cell_type in cell_types
            }
            max_cell, max_mean = max(means.items(), key=lambda x: x[1])
            enrichment = max_mean / (sum(means.values()) - max_mean)

            enriched[species][row[gene_col_name]] = max_cell, enrichment

    try:
        with open(cache.ENRICHMENT_CACHE, "wb") as f:
            pickle.dump(enriched, f)
    except:
        pass

    return enriched


def get_enrichments(species, add_mappings=True, cutoff=2.5, **kwargs):
    enrichments = {
        key: val
        for key, val in build_enrichment_table(**kwargs)[species].items()
        if val[1] >= cutoff
    }

    for index, row in cache.get_mapping_data().iterrows():
        if index not in enrichments:
            continue

        for syn in row["Synonyms"].split("|"):
            if syn in enrichments or syn == "-":
                continue

            enrichments[syn] = enrichments[index]

    return enrichments