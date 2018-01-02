
import logging
import pickle

import numpy as np

from . import cache, mapping

# Taken from
LOGGER = logging.getLogger("brainrnaseq")


CELL_TYPE_COLS = {
    "Human": {
        "Astrocyte": [
            '8yo',
            '13yo', '16yo', '21yo.1', '22yo.1', '35yo', '47yo', '51yo', '53yo',
            '60yo', '63yo - 1', '63yo - 2',
        ],
        "Neuron": [
            "25yo",
        ],
        "OPC": [
            '22yoGC', '63yoGC - 1',
            '63yo GC - 2', '47yoO4', '63yoO4',
        ],
        "New Oligodendrocytes": [
            '22yoGC', '63yoGC - 1',
            '63yo GC - 2', '47yoO4', '63yoO4',
        ],
        "Myelinating Oligodendrocytes": [
            '22yoGC', '63yoGC - 1',
            '63yo GC - 2', '47yoO4', '63yoO4',
        ],
        "Microglia": [
            '45yo', '51yo.1', '63yo',
        ],
        "Endothelia": [
            "13yo.1",
            "47yo.1",
        ],
    },
    "Mouse": {
        "Astrocyte": [
            "FACS - p69",
            "FACS p70",
            "1 month",
            "4 months",
            "7 months",
            "9 months",
        ],
        "Neuron": [
            "Neuron 3",
            "Neuron 4",
        ],
        "OPC": [
            "Oligodendrocyte precursor cell 3",
            "Oligodendrocyte precursor cell 4",
        ],
        "New Oligodendrocytes": [
            "Newly formed oligodendrocyte 3",
            "Newly formed oligodendrocyte 4",
        ],
        "Myelinating Oligodendrocytes": [
            "Myelinating oligodendrocyte 4",
            "Myelinating oligodenrocyte 5",
        ],
        "Microglia": [
            "Microglia 1",
            "Microglia 2",
        ],
        "Endothelia": [
            "Endo 1",
            "Endo 2",
        ],
    },
}

CELL_TYPES = [
    "Astrocyte",
    "Neuron",
    "OPC",
    "New Oligodendrocytes",
    "Myelinating Oligodendrocytes",
    "Microglia",
    "Endothelia",
]
DEFAULT_CELL_TYPES = CELL_TYPES[:2] + CELL_TYPES[4:]


def build_enrichment_table(cell_types=None, force=False):
    if not force:
        try:
            with open(cache.ENRICHMENT_CACHE, "rb") as f:
                return pickle.load(f)
        except:
            pass

    if cell_types is None:
        cell_types = DEFAULT_CELL_TYPES

    cache.get_barres_seq_data(force=force)

    LOGGER.info("Calculating cell type enrichment")
    enriched = {}

    for species, data in SPECIES_DATA.items():
        enriched[species] = {}

        for col in data.columns:
            if hasattr(col, "lower") and col.lower() == "gene":
                gene_col_name = col
                break

        data = data.drop_duplicates(subset=gene_col_name)

        for _, row in data.iterrows():
            if any(
                isinstance(row[col], str)
                for cols in CELL_TYPE_COLS[species].values()
                for col in cols
            ):
                continue

            means = {
                cell_type: np.nanmean(np.array(row[cols].values, dtype=float))
                for cell_type, cols in CELL_TYPE_COLS[species].items()
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
