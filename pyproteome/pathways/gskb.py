
import logging
import requests

import pandas as pd

import pyproteome as pyp

LOGGER = logging.getLogger('pyproteome.gskb')

GSKB_URL = (
    'http://ge-lab.org/gskb/2-MousePath/mGSKB_Entrez.gmt'
)


@pyp.utils.memoize
def get_gskb_pathways(species):
    '''
    Download gene sets from GSKB.

    Parameters
    ----------
    species : str

    Returns
    -------
    df : :class:`pandas.DataFrame`, optional
    '''
    LOGGER.info('Fetching GSKB pathways')

    url = GSKB_URL
    r = requests.get(url, stream=True)
    r.raise_for_status()

    def _get_data(line):
        line = line.decode('windows-1252')
        name, _, genes = line.split('\t', 2)

        genes = set(i for i in genes.split('\t') if i)

        return name, genes

    pathways_df = pd.DataFrame(
        data=[
            _get_data(line)
            for ind, line in enumerate(r.iter_lines())
            if ind > 0
        ],
        columns=['name', 'set'],
    )

    return pathways_df
