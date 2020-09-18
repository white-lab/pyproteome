
import gzip
import io
import logging
import re
import requests

import pandas as pd

import pyproteome as pyp
import brainrnaseq as brs

LOGGER = logging.getLogger('pyproteome.pathwayscommon')

PATHWAYS_COMMON_URL = (
    'http://www.pathwaycommons.org/archives/PC2/v9/'
    'PathwayCommons9.All.hgnc.gmt.gz'
)


@pyp.utils.memoize
def get_pathway_common(species):
    '''
    Download gene sets from Pathway Commons.

    Parameters
    ----------
    species : str

    Returns
    -------
    df : :class:`pandas.DataFrame`, optional
    '''
    LOGGER.info('Fetching Pathways Common')

    url = PATHWAYS_COMMON_URL
    r = requests.get(url, stream=True)
    r.raise_for_status()

    name_re = re.compile(
        'name: (.+); datasource: (.+); organism: (.+); idtype: (.+)'
    )

    def _get_data(line):
        line = line.decode('utf-8')
        _, name, genes = line.split('\t', 2)
        name = name_re.match(name)

        name = {
            'name': name.group(1),
            'datasource': name.group(2),
            'organism': name.group(3),
            'id_type': name.group(4),
        }

        assert int(name['organism']) == 9606

        genes = set(
            brs.mapping.get_entrez_mapping(i, species=species)
            for i in genes.split('\t')
        )

        return name['name'], genes

    pathways_df = pd.DataFrame(
        data=[
            _get_data(line)
            for line in gzip.GzipFile(fileobj=io.BytesIO(r.content))
        ],
        columns=['name', 'set'],
    )

    return pathways_df
