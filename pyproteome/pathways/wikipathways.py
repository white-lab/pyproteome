
import io
import logging
import requests
import re
import zipfile
import xml.etree.ElementTree as ET

import pandas as pd

import pyproteome as pyp

WIKIPATHWAYS_GMT_URL = (
    'http://data.wikipathways.org/{date}/gmt/'
    'wikipathways-{date}-gmt-{species}.gmt'
)
WIKIPATHWAYS_GMPL_URL = (
    'http://data.wikipathways.org/{date}/gpml/'
    'wikipathways-{date}-gpml-{species}.zip'
)
WIKIPATHWAYS_NS = {
    'wp': 'http://pathvisio.org/GPML/2013a',
}
RE_WIKIPATHWAYS = re.compile(
    r'parent=([^;]+); position=([^;]+); ptm=([^;]+); direction=([^;]+)'
)

LOGGER = logging.getLogger('pyproteome.wikipathways')


@pyp.utils.memoize
def _get_wp_date():
    return '20180710'


@pyp.utils.memoize
def get_wikipathways(species):
    '''
    Download gene sets from WikiPathways.

    Parameters
    ----------
    species : str

    Returns
    -------
    df : :class:`pandas.DataFrame`, optional
    '''
    LOGGER.info('Fetching WikiPathways')

    species = pyp.species.INV_ORGANISM_MAPPING.get(species, species)

    url = WIKIPATHWAYS_GMT_URL.format(
        date=_get_wp_date(),
        species='_'.join(species.split(' ')),
    )
    response = requests.get(url, stream=True)
    response.raise_for_status()

    def _get_data(line):
        line = line.decode('utf-8')
        name, _, genes = line.split('\t', 2)
        name, _, _, spec = name.split('%')
        assert species == spec
        return name, set(i for i in genes.split('\t'))

    pathways_df = pd.DataFrame(
        data=[
            _get_data(line)
            for line in response.iter_lines()
        ],
        columns=['name', 'set'],
    )

    return pathways_df


@pyp.utils.memoize
def get_wikipathways_psites(species):
    '''
    Download phospho sets from WikiPathways.

    Parameters
    ----------
    species : str

    Returns
    -------
    df : :class:`pandas.DataFrame`, optional
    '''
    LOGGER.info('Fetching WikiPathways')

    species = pyp.species.INV_ORGANISM_MAPPING.get(species, species)

    url = WIKIPATHWAYS_GMPL_URL.format(
        date=_get_wp_date(),
        species='_'.join(species.split(' ')),
    )
    response = requests.get(url, stream=False)
    response.raise_for_status()

    z = zipfile.ZipFile(io.BytesIO(response.content))
    LOGGER.info('Parsing WikiPathways phosphosites')

    def _process_site(sname):
        sname = re.sub('ser', 'S', sname)
        sname = re.sub('thr', 'T', sname)
        sname = re.sub('tyr', 'Y', sname)
        return sname

    def _get_data(name):
        f = z.open(name)
        root = ET.fromstring(f.read())
        sites = root.findall(
            'wp:State/wp:Comment',
            WIKIPATHWAYS_NS,
        )
        # print(sites)
        matches = [
            RE_WIKIPATHWAYS.match(i.text)
            for i in sites
            if i.text
        ]
        matches = [
            (
                m.group(1),
                _process_site(m.group(2)),
                m.group(4),
            )
            for m in matches
            if m and m.group(3) == 'p'
        ]
        up = [i for i in matches if i[2] == 'u']
        down = [i for i in matches if i[2] == 'd']
        # parent=Q9NQB0; position=thr212; ptm=p; direction=d

        tmp = name.replace('_', ' ').rsplit('.', 1)[0].split(' ', 1)[1]
        title, id = tmp.split(' WP')
        id = 'WP' + id

        return (
            title,
            id,
            set(
                '{},{}-p'.format(i[0], i[1])
                for i in up
            ),
            set(
                '{},{}-p'.format(i[0], i[1])
                for i in down
            ),
        )

    data = [
        _get_data(name)
        for name in z.namelist()
    ]
    data = [
        i
        for i in data
        if i[1] or i[2]
    ]

    pathways_df = pd.DataFrame(
        data=data,
        columns=['name', 'WikiPathways ID', 'upregulated', 'downregulated'],
    )
    pathways_df = pathways_df[
        pathways_df['upregulated'].apply(len) +
        pathways_df['downregulated'].apply(len) > 5
    ]

    return pathways_df
