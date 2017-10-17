
import logging
import os
import pickle
import requests

import numpy as np
import pandas as pd

# Taken from
LOGGER = logging.getLogger("brainrnaseq")

BARRES_SEQ_URL = (
    "https://web.stanford.edu/group/barres_lab/brainseq2/"
    "TableS4-HumanMouseMasterFPKMList.xlsx"
)
DIR = os.path.abspath(os.path.split(__file__)[0])
BARRES_DATA_NAME = "TableS4-HumanMouseMasterFPKMList.xlsx"
CACHE_DIR = os.path.join(DIR, "cache")
BARRES_SEQ_PATH = os.path.join(CACHE_DIR, BARRES_DATA_NAME)
ENRICHMENT_CACHE = os.path.join(CACHE_DIR, "enrichment_cache.pickle")
SPECIES_DATA = {}

try:
    os.makedirs(CACHE_DIR)
except:
    pass

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

GENE_MAPPING = {
    "Champ1": "Zfp828",
    "Dcdc2":  "Dcdc2a",
    "Map":    "Mtap",
    "Map1a":  "Mtap1a",
    "Map1b":  "Mtap1b",
    "Map1s":  "Mtap1s",
    "Map2":   "Mtap2",
    "Map4":   "Mtap4",
    "Map6":   "Mtap6",
    "Map7":   "Mtap7",
    "Map7d1": "Mtap7d1",
    "Map7d2": "Mtap7d2",
    "Map7d3": "Mtap7d3",
    "Map9":   "Mtap9",
    "Prkcg":  "Prkcc",
    "Znf148": "Zfp148",
    "Dmtn":   "Epb4.9",
    "Skt":    "Etl4",
    "Plppr4": "Ier3",
    "Plpp4":  "Ppapdc1a",
    "Lgalsl": "1110067D22Rik",
    "Sptbn1": "Spnb2",
    "Kiaa1109": "4932438A13Rik",
    "Carmil2": "Rltpr",
    "Apmap":   "2310001A20Rik",
    "Lurap1":  "1520402A15Rik",
    "Tp53bp1": "Trp53bp1",
    "Cep128":  "4930534B04Rik",
    "Mb21d2":  "1600021P15Rik",
    "Cep41":   "Tsga14",
    "Plekhd1": "3830431G21Rik",
    "Camsap2": "Camsap1l1",
    "Lamtor1": "2400001E08Rik",
    "Gpi":     "Gpi1",
    "Fer":     "Fert2",
    "Nyap2":   "9430031J16Rik",
    "Ctc1":    "1500010J02Rik",
    "Kiaa2022": "C77370",
    "Tubb4a":   "Tubb4",
    "Pea15":    "Pea15a",
    "Fam208b":  "BC016423",
    "Lrif1":    "4933421E11Rik",
    "Dnal1":    "Dnalc1",
    "Znf687":   "Zfp687",
    "Smg9":     "1500002O20Rik",
}


def get_barres_seq_data(force=False):
    global SPECIES_DATA

    if force or not os.path.exists(BARRES_SEQ_PATH):
        LOGGER.info("Downloading Barres RNA Seq Data")
        response = requests.get(BARRES_SEQ_URL, stream=True)
        response.raise_for_status()

        with open(BARRES_SEQ_PATH, mode="wb") as f:
            for block in response.iter_content(1024):
                f.write(block)

    LOGGER.info("Reading Barres RNA Seq Data")
    SPECIES_DATA = {
        "Human": pd.read_excel(
            BARRES_SEQ_PATH, sheetname="Human data only",
            skiprows=[0],
        ),
        "Mouse": pd.read_excel(
            BARRES_SEQ_PATH, sheetname="Mouse data only",
            skiprows=[0],
        ),
    }


def build_enrichment_table(cell_types=None, force=False):
    if not force:
        try:
            with open(ENRICHMENT_CACHE, "rb") as f:
                return pickle.load(f)
        except:
            pass

    if cell_types is None:
        cell_types = DEFAULT_CELL_TYPES

    get_barres_seq_data(force=force)

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
        with open(ENRICHMENT_CACHE, "wb") as f:
            pickle.dump(enriched, f)
    except:
        pass

    return enriched


def get_enrichments(species, cutoff=2.5, **kwargs):
    return {
        key: val
        for key, val in build_enrichment_table(**kwargs)[species].items()
        if val[1] >= cutoff
    }
