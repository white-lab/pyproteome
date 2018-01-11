
import logging
import io
import os
import pickle
import requests
import zlib

import pandas as pd

# Taken from
LOGGER = logging.getLogger("brainrnaseq.cache")
DIR = os.path.abspath(os.path.split(__file__)[0])
CACHE_DIR = os.path.join(DIR, "cache")

# BARRES_SEQ_URL = (
#     "https://web.stanford.edu/group/barres_lab/brainseq2/"
#     "TableS4-HumanMouseMasterFPKMList.xlsx"
# )
BARRES_SEQ_URL = (
    "https://github.com/white-lab/pyproteome-data/blob/master/brainrnaseq/"
    "TableS4-HumanMouseMasterFPKMList.xlsx?raw=true"
)
BARRES_DATA_NAME = "TableS4-HumanMouseMasterFPKMList.xlsx"

BARRES_SEQ_PATH = os.path.join(CACHE_DIR, BARRES_DATA_NAME)

# MAPPING_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/"
# HUMAN_MAPPING_URL = MAPPING_URL + "Homo_sapiens.gene_info.gz"
# MOUSE_MAPPING_URL = MAPPING_URL + "Mus_musculus.gene_info.gz"
MAPPING_URL = (
    "https://github.com/white-lab/pyproteome-data/raw/master/brainrnaseq/"
)
HUMAN_MAPPING_URL = MAPPING_URL + "Homo_sapiens.gene_info.gz"
MOUSE_MAPPING_URL = MAPPING_URL + "Mus_musculus.gene_info.gz"

MAPPING_URLS = {
    "Human": HUMAN_MAPPING_URL,
    "Mouse": MOUSE_MAPPING_URL
}

ENRICHMENT_CACHE = os.path.join(CACHE_DIR, "enrichment_cache.pickle")
MAPPING_CACHE = os.path.join(CACHE_DIR, "mapping_cache.pickle")
SPECIES_DATA = {}
MAPPING_DATA = None


try:
    os.makedirs(CACHE_DIR)
except:
    pass


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
            BARRES_SEQ_PATH,
            sheet_name="Human data only",
            skiprows=[0],
        ).iloc[1:],
        "Mouse": pd.read_excel(
            BARRES_SEQ_PATH,
            sheet_name="Mouse data only",
            skiprows=[0],
        ),
    }


def fetch_mapping_data(url):
    LOGGER.info("Fetching mapping data from {}".format(url))

    response = requests.get(url, stream=True)
    response.raise_for_status()

    dec = zlib.decompressobj(32 + zlib.MAX_WBITS)
    csv = b""

    for block in response.iter_content(1024):
        csv += dec.decompress(block)

    df = pd.read_csv(
        io.StringIO(csv.decode("utf-8")),
        sep="\t",
    )

    LOGGER.info("Read mapping info for {} genes".format(df.shape[0]))

    df = df.set_index("Symbol")

    return df


def get_mapping_data(species="Mouse", force=False):
    global MAPPING_DATA
    global MAPPING_CACHE

    if MAPPING_DATA is None:
        if not force:
            try:
                with open(MAPPING_CACHE, "rb") as f:
                    MAPPING_DATA = pickle.load(f)
            except:
                MAPPING_DATA = {}
        else:
            MAPPING_DATA = {}

    assert species in MAPPING_URLS

    if species in MAPPING_DATA:
        return MAPPING_DATA[species]

    MAPPING_DATA[species] = fetch_mapping_data(
        MAPPING_URLS[species]
    )

    try:
        with open(MAPPING_CACHE, "wb") as f:
            pickle.dump(MAPPING_DATA, f)
    except Exception as e:
        LOGGER.warning(
            "Unable to save mapping information to cache: {}"
            .format(e)
        )

    return MAPPING_DATA[species]
