
from collections import OrderedDict
import io
import logging
import requests
import os
import tarfile
import uuid

import pyproteome as pyp
import brainrnaseq as brs
import pandas as pd
import numpy as np

LOGGER = logging.getLogger('pathways.photon_ptm')

try:
    from genemap.mappers import EnsemblMapper
except ImportError:
    pass


parameters = {
    'activity': {
        'min_size': 4,
        'permutations': 1000,
        'side': 'greater',
    },
    'anat': {
        'alpha': 0.25,
        'anchor': -1,
    },
    'go': {
        'max_category_size': 500,
        'min_category_size': 5,
        'max_category_depth': 10,
    },
    'ppi-network': {
        'confidence': 0.5,
        'degree_threshold': 150,
    }
}


@pyp.utils.memoize
def _get_anat(dir, force=False):
    if os.path.exists(os.path.join(dir, 'db')) and not force:
        return

    url = 'http://cs.tau.ac.il/~jandanielr/db.tar.gz'

    LOGGER.info('Fetching PHOTON database from {} to {}'.format(url, dir))

    r = requests.get(url, stream=True)
    r.raw.decode_content = True
    r.raise_for_status()

    tar_file = tarfile.open(
        fileobj=io.BytesIO(r.content),
        mode='r',
    )
    tar_file.extractall(dir)

    return tar_file


@pyp.utils.memoize
def _map_gene(mapper, symbol_mapping, gene, species):
    symbol = gene
    entrez = brs.mapping.get_entrez_mapping(symbol, species=species)

    if species.lower().replace(' ', '_') not in ['homo_sapiens', 'human']:
        entrez = [i for i in mapper.map_ids([entrez]) if i]

        if not entrez:
            return None, None

        entrez = int(entrez[0])
        symbol = symbol_mapping.loc[entrez]

    return entrez, symbol


def _get_templates(template_dir):
    pyp.utils.makedirs(template_dir)

    url = (
        'https://raw.githubusercontent.com/jdrudolph/photon/'
        'c73d1eb7f5e7cab86031e056350c7b09fb5e1d51/templates/result.html'
    )
    r = requests.get(url)
    r.raise_for_status()

    with open(os.path.join(template_dir, 'result.html'), 'wb') as f:
        f.write(r.content)


def photon(ds, folder_name=None, write_output=False):
    '''
    Run PHOTON algorithm on a data set to find functional phosphorylation sites
    using protein-protein interaction networks.

    Parameters
    ----------
    ds : :class:`pyproteome.data_sets.DataSet`
    folder_name : str, optional

    Returns
    -------
    out_path : str
        Path to results directory.
    '''
    import phos
    import phos.defaults
    import phos.pipeline

    ds = ds.filter(fn=lambda x: len(x['Proteins']) < 2)
    ds = ds.filter(fn=lambda x: len(x['Scan']) >= 2)
    ds = ds.filter(mod='Phospho')

    species = list(ds.species)[0]

    from_name = '{}{}'.format(
        species.split(' ')[0][0],
        species.split(' ')[1],
    ).lower()

    species = species.replace(' ', '_')

    mapper = EnsemblMapper(
        from_type='entrez',
        to_type='entrez',
        from_organism=from_name,
        to_organism='hsapiens',
    )

    symbol_mapping = brs.cache.get_mapping_data(species='Homo sapiens')
    symbol_mapping['Symbol'] = symbol_mapping.index
    symbol_mapping = symbol_mapping.set_index('GeneID')['Symbol']

    def _get_phos_data(psms):
        hit = 0
        miss = 0

        for _, row in psms.iterrows():
            gene = row['Proteins'].genes[0]

            entrez, symbol = _map_gene(mapper, symbol_mapping, gene, species)

            if not entrez:
                # print(gene, entrez, symbol)
                miss += 1
                continue

            hit += 1

            for mod in row['Modifications'].get_mods('Phospho'):
                yield pd.Series(OrderedDict([
                    ('GeneID', entrez),
                    ('Amino.Acid', mod.letter),
                    ('Position', 1 + mod.abs_pos[0]),
                    ('avg', np.log2(row['Fold Change'])),
                    ('Symbol', symbol),
                ]))

        print(hit, miss)

    LOGGER.info('Generated data frame: {}, {}'.format(ds.name, ds.shape))

    df = pd.DataFrame(
        list(_get_phos_data(ds.psms))
    ).dropna()

    LOGGER.info('Generated data frame: {}'.format(df.shape))

    df = df.sort_values('avg', ascending=False)

    name = str(uuid.uuid4())

    _parameters = parameters.copy()
    _parameters['anat']['anchor'] = 1950

    # XXX: Better directory for files that Photon must download for its use?
    dir = os.path.join(pyp.paths.FIGURES_DIR)
    defaults = phos.defaults.make_defaults(dir)
    _get_anat(dir)

    template_dir = os.path.join(defaults['root'], 'templates')
    _get_templates(template_dir)

    folder_name = pyp.utils.make_folder(
        data=ds,
        folder_name=folder_name,
        sub=os.path.join('Photon', name),
    )

    csv_path = os.path.join(folder_name, 'results.csv')

    with open(csv_path, 'w') as csv_file:
        df.to_csv(csv_file, index=False)

    if write_output:
        phos.pipeline.run(
            name,
            csv_path,
            _parameters,
            template_dir,
            folder_name,
            defaults['db'],
        )

        LOGGER.info('Wrote results to: {}'.format(folder_name))

        return folder_name
    else:
        exp, scores, subnet, go_scores, predictions = phos.pipeline._run(
            name,
            csv_path,
            _parameters,
            # template_dir,
            # folder_name,
            defaults['db'],
        )

        return exp, scores, subnet, go_scores, predictions
