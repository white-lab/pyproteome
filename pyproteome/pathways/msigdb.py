
import logging
import requests

import pandas as pd

import pyproteome as pyp


MSIGDB_URL = (
    'https://raw.githubusercontent.com/white-lab/pyproteome-data'
    '/master/msigdb/'
)
MSIGDB_VERSION = 'v7.1'
MSIGDB_FILES = tuple([
    i.format(MSIGDB_VERSION)
    for i in [
        'h.all.{}.entrez.gmt',
        # 'c1.all.{}.entrez.gmt',
        # 'c2.all.{}.entrez.gmt',
        # 'c2.cgp.{}.entrez.gmt',
        # 'c2.cp.biocarta.{}.entrez.gmt',
        # 'c2.cp.kegg.{}.entrez.gmt',
        # 'c2.cp.reactome.{}.entrez.gmt',
        # 'c2.cp.{}.entrez.gmt',
        # 'c3.all.{}.entrez.gmt',
        # 'c3.mir.{}.entrez.gmt',
        # 'c3.tft.{}.entrez.gmt',
        # 'c4.all.{}.entrez.gmt',
        # 'c4.cgn.{}.entrez.gmt',
        # 'c4.cm.{}.entrez.gmt',
        # 'c5.all.{}.entrez.gmt',
        # 'c5.bp.{}.entrez.gmt',
        # 'c5.cc.{}.entrez.gmt',
        # 'c5.mf.{}.entrez.gmt',
        # 'c6.all.{}.entrez.gmt',
        # 'c7.all.{}.entrez.gmt',
    ]
])

LOGGER = logging.getLogger('pyproteome.msigdb')

try:
    from genemap.mappers import EnsemblMapper
except ImportError:
    pass


@pyp.utils.memoize
def get_msigdb_pathways(species, remap=None):
    '''
    Download gene sets from MSigDB. Currently downloads v7.0 of the gene
    signature repositories.

    Parameters
    ----------
    species : str

    Returns
    -------
    df : :class:`pandas.DataFrame`, optional
    '''
    LOGGER.info('Fetching MSigDB pathways')

    def _get_requests():
        for file in MSIGDB_FILES:
            url = MSIGDB_URL + file

            LOGGER.info('Fetching {}'.format(url))

            response = requests.get(url, stream=True)
            response.raise_for_status()

            yield response

    def _get_data(line):
        line = line.decode('utf-8')
        name, _, genes = line.split('\t', 2)
        # name, _, _, spec = name.split('%')
        # assert species == spec
        return name, set(i for i in genes.split('\t'))

    pathways_df = pd.DataFrame(
        data=[
            _get_data(line)
            for response in _get_requests()
            for line in response.iter_lines()
        ],
        columns=['name', 'set'],
    )

    if remap and species not in ['Homo sapiens']:
        to_name = '{}{}'.format(
            species.split(' ')[0][0],
            species.split(' ')[1],
        ).lower()

        LOGGER.info('Remapping MSigDB to {} ({})'.format(species, to_name))

        mapper = EnsemblMapper(
            from_type='entrez',
            to_type='entrez',
            from_organism='hsapiens',
            to_organism=to_name,
        )
        pathways_df['set'] = pathways_df['set'].apply(
            lambda row: set(mapper.map_ids(row))
        )

    return pathways_df
