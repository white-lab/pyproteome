
import gzip
import io
import logging
import requests

import pandas as pd

import pyproteome as pyp

LOGGER = logging.getLogger('pyproteome.phosphosite')

PSP_REGULATORY_URL = (
    'https://www.phosphosite.org/downloads/Regulatory_sites.gz'
)
PSP_SITE_MAPPING_URL = (
    'https://www.phosphosite.org/downloads/Phosphorylation_site_dataset.gz'
)


@pyp.utils.memoize
def get_phosphomap_data():
    '''
    Fetch mapping between phosphorylation sites of different species.

    Returns
    -------
    df : :class:`pandas.DataFrame`
    '''
    LOGGER.info('Fetching Phosphosite Plus mapping data')

    url = PSP_SITE_MAPPING_URL

    r = requests.get(url, stream=True)
    r.raise_for_status()

    gz = gzip.GzipFile(fileobj=io.BytesIO(r.content))

    return pd.read_csv(gz, skiprows=[0, 1, 2], sep='\t')


@pyp.utils.memoize
def get_phosphoreg_data():
    '''
    Fetch Phosphosite Plus regulation data.

    Returns
    -------
    df : :class:`pandas.DataFrame`
    '''
    LOGGER.info('Fetching Phosphosite Plus regulation data')

    url = PSP_REGULATORY_URL

    r = requests.get(url, stream=True)
    r.raise_for_status()

    gz = gzip.GzipFile(fileobj=io.BytesIO(r.content))

    return pd.read_table(gz, skiprows=[0, 1, 2], sep='\t', usecols=range(21))


@pyp.utils.memoize
def get_phosphosite(species, remap=False):
    '''
    Download phospho sets from PhophoSite Plus.

    Parameters
    ----------
    species : str
    remap : bool, optional

    Returns
    -------
    df : :class:`pandas.DataFrame`, optional
    '''
    LOGGER.info('Getting phosphosite data for {}'.format(species))

    species = pyp.species.ORGANISM_MAPPING.get(species, species)

    psp_data = pyp.motifs.phosphosite.get_data()

    if remap:
        psp_data = _remap_psp(
            psp_data, species,
            acc_col='SUB_ACC_ID',
            mod_col='SUB_MOD_RSD',
            org_col='SUB_ORGANISM',
        )

    psp_data = psp_data[psp_data['SUB_ORGANISM'] == species]

    return pd.DataFrame(
        [
            (
                kinase,
                set(
                    psp_data[
                        psp_data['KINASE'] == kinase
                    ].apply(
                        lambda x:
                        ','.join([
                            x['SUB_ACC_ID'].split('-')[0]
                            if isinstance(x['SUB_ACC_ID'], str) else
                            '',
                            x['SUB_MOD_RSD'] + ('' if remap else '-p'),
                        ]),
                        axis=1,
                    )
                ),
                set(),
            )
            for kinase in set(psp_data['KINASE'])
        ],
        columns=['name', 'up_set', 'down_set']
    )


@pyp.utils.memoize
def get_phosphosite_regulation(species, remap=False):
    '''
    Download phospho sets from PhophoSite Plus.

    Parameters
    ----------
    species : str
    remap : bool, optional

    Returns
    -------
    df : :class:`pandas.DataFrame`, optional
    '''
    LOGGER.info('Getting phosphosite regulation data for {}'.format(species))

    species = pyp.species.ORGANISM_MAPPING.get(species, species).lower()

    psp_data = get_phosphoreg_data()

    if remap:
        psp_data = _remap_psp(
            psp_data, species,
            set_col='ON_PROCESS',
            acc_col='ACC_ID',
            mod_col='MOD_RSD',
            org_col='ORGANISM',
            append_mod='',
        )

    psp_data = psp_data[psp_data['ORGANISM'] == species]

    paths = set(
        proc.strip()
        for row in psp_data['ON_PROCESS']
        if pd.notna(row)
        for proc in row.split(';')
        if proc.strip()
    )

    return pd.DataFrame(
        [
            (
                path,
                set(
                    psp_data[
                        psp_data['ON_PROCESS'].apply(
                            lambda x:
                            pd.notna(x) and any([
                                path in i.strip()
                                for i in x.split(';')
                            ])
                        )
                    ].apply(
                        lambda x:
                        ','.join([
                            x['ACC_ID'].split('-')[0],
                            x['MOD_RSD'],
                        ]),
                        axis=1,
                    )
                ),
                set(),
            )
            for path in paths
        ],
        columns=['name', 'up_set', 'down_set']
    )


def _remap_psp(
    psp, species,
    set_col='KINASE',
    acc_col='ACC_ID',
    mod_col='MOD_RSD',
    org_col='ORGANISM',
    append_mod='-p',
    mapping=None,
):
    LOGGER.info('Remapping sites to species: {}'.format(species))

    if mapping is None:
        mapping = get_phosphomap_data()
    
    mapping = mapping[['ACC_ID', 'MOD_RSD', 'ORGANISM', 'SITE_GRP_ID']]

    mod_mapping = mapping[
        mapping['ORGANISM'] != species
    ].set_index(
        ['ACC_ID', 'MOD_RSD', 'ORGANISM']
    ).sort_index()
    site_mapping = mapping[
        mapping['ORGANISM'] == species
    ].set_index(
        'SITE_GRP_ID'
    ).sort_index()
    del mapping

    new_index = [org_col, set_col, acc_col, mod_col]

    def _remap(row):
        kinase, acc, mod, old_species = row[
            [set_col, acc_col, mod_col, org_col]
        ]
        mod += append_mod

        if old_species != species:
            # Remap the phosphorylation site if possible
            try:
                site = mod_mapping.loc(axis=0)[acc, mod, old_species]
            except KeyError:
                pass
            else:
                site = site['SITE_GRP_ID']

                if hasattr(site, 'iloc'):
                    site = site.iloc[0]

                try:
                    re_map = site_mapping.loc[site]
                except KeyError:
                    pass
                else:
                    if len(re_map.shape) > 1:
                        re_map = re_map.iloc[0]

                    acc, mod = re_map[['ACC_ID', 'MOD_RSD']]
                    old_species = species

        return pd.Series([
            old_species,
            kinase,
            acc,
            mod,
        ], index=new_index)

    return psp.apply(_remap, axis=1)
