'''
This file includes functions for downloading kinase-substrate associations from
PhosphoSite Plus (https://www.phosphosite.org/).
'''

from io import BytesIO
import gzip
import os
import requests

import numpy as np
import pandas as pd

import pyproteome as pyp
from . import motif, logo


DATA_URL = 'https://www.phosphosite.org/downloads/Kinase_Substrate_Dataset.gz'


@pyp.utils.memoize
def get_data():
    '''
    Download the Kinase-Substrate Dataset from Phosphosite Plus.

    Returns
    -------
    df : :class:`pandas.DataFrame`
    '''
    data = requests.get(DATA_URL, stream=True)
    content = BytesIO(data.content)

    with gzip.GzipFile(fileobj=content) as f:
        df = pd.read_csv(f, skiprows=range(2), sep='\t')

    return df


def generate_logos(
    species,
    kinases=None,
    min_foreground=10,
):
    '''
    Generate logos for all kinases documented on Phosphosite Plus.

    Parameters
    ----------
    species : str
        Species name (i.e. 'Human' or 'Homo sapiens')
    kinases : list of str, optional
    min_foreground : int, optional
        Minimum number of substrates needed for logo generation.
    '''
    species = pyp.species.ORGANISM_MAPPING.get(species, species).lower()

    df = get_data()
    df = df[
        np.logical_and(
            df['KIN_ORGANISM'].apply(lambda x: x.lower()) == species,
            df['SUB_ORGANISM'].apply(lambda x: x.lower()) == species,
        )
    ]

    if kinases is None:
        kinases = sorted(set(df['KINASE']))

    figs = []

    for kinase in kinases:
        fore = list(df[df['KINASE'] == kinase]['SITE_+/-7_AA'])

        if len(fore) < min_foreground:
            continue

        f = logo.logo(
            fore=fore,
            back=list(df['SITE_+/-7_AA']),
            title=kinase,
        )[0]

        figs.append(f)

    return figs


def enriched(data, species=None):
    df = get_data()

    if species:
        df = df[
            np.logical_and(
                df['KIN_ORGANISM'] == species,
                df['SUB_ORGANISM'] == species,
            )
        ]

    return df[
        df['SITE_+/-7_AA'].isin(
            motif.generate_n_mers(
                data['Sequence'],
                fill_left='_',
                fill_right='_',
                mods=[(None, 'Phospho')],
            )
        )
    ].style.set_table_styles(
        [
            {'selector': 'th:first-child', 'props': [('display', 'none')]},
        ]
    )
