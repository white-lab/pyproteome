
from io import BytesIO
import gzip
import requests

import numpy as np
import pandas as pd

from . import motif

DATA_URL = "https://www.phosphosite.org/downloads/Kinase_Substrate_Dataset.gz"
DATA_CACHE = None


def get_data():
    global DATA_CACHE

    if DATA_CACHE is not None:
        return DATA_CACHE

    # r = requests.Session()

    data = requests.get(DATA_URL, stream=True)
    content = BytesIO(data.content)

    with gzip.GzipFile(fileobj=content) as f:
        df = pd.read_csv(f, skiprows=range(2), sep="\t")

    DATA_CACHE = df

    return df


def enriched(data, species=None):
    df = get_data()

    if species:
        df = df[
            np.logical_and(
                df["KIN_ORGANISM"] == species,
                df["SUB_ORGANISM"] == species,
            )
        ]

    return df[
        df["SITE_+/-7_AA"].isin(
            motif.generate_n_mers(
                data["Sequence"],
                fill_left="_",
                fill_right="_",
                letter_mod_types=[(None, "Phospho")],
            )
        )
    ]
