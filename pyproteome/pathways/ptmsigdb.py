
import os
import requests

import pandas as pd

import pyproteome as pyp


PTMSIGDB_URL = (
    'https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/'
    'ptmsigdb/ptm.sig.db.all.uniprot.{name}.v1.9.0.gmt'
)


@pyp.utils.memoize
def get_ptmsigdb(species):
    '''
    Download phospho sets for PTMSigDB.

    Parameters
    ----------
    species : str

    Returns
    -------
    df : :class:`pandas.DataFrame`
    '''
    species = pyp.species.ORGANISM_MAPPING.get(species, species).lower()

    assert species in ['human', 'mouse', 'rat']

    def _get_data(line):
        line = line.split('\t')
        sites = set(
            i.strip().replace(';', ',')
            for i in line[2:]
        )
        up_sites = set(
            i.rsplit(',', 1)[0]
            for i in sites
            if i.endswith('u')
        )
        down_sites = set(
            i.rsplit(',', 1)[0]
            for i in sites
            if i.endswith('d')
        )
        title, description = line[:2]
        return (
            title,
            description,
            up_sites,
            down_sites
        )

    url = PTMSIGDB_URL.format(name=species.lower())

    response = requests.get(url, stream=True)
    response.raise_for_status()

    data = [
        _get_data(line)
        for line in response.content.decode('utf-8').split('\n')
        if line.strip()
    ]

    pathways_df = pd.DataFrame(
        data=data,
        columns=['name', 'Description', 'up_set', 'down_set'],
    )

    return pathways_df[['name', 'up_set', 'down_set']]
