
from collections import OrderedDict
import io
import logging
import requests
import os
import tarfile
import tempfile
import uuid

import pyproteome as pyp
import brainrnaseq as brs
import pandas as pd

LOGGER = logging.getLogger("pathways.photon_ptm")

try:
    from genemap.mappers import EnsemblMapper
except ImportError:
    pass


parameters = {
    "activity": {
        "min_size": 4,
        "permutations": 1000,
        "side": "greater",
    },
    "anat": {
        "alpha": 0.25,
        "anchor": -1,
    },
    "go": {
        "max_category_size": 500,
        "min_category_size": 5,
        "max_category_depth": 10,
    },
    "ppi-network": {
        "confidence": 0.5,
        "degree_threshold": 150,
    }
}


@pyp.utils.memoize
def _get_anat(dir, force=False):
    if os.path.exists(os.path.join(dir, "db")) and not force:
        return

    url = "http://cs.tau.ac.il/~jandanielr/db.tar.gz"

    LOGGER.info("Fetching PHOTON database from {} to {}".format(url, dir))

    r = requests.get(url, stream=True)
    r.raw.decode_content = True
    r.raise_for_status()

    # gzip_file = gzip.GzipFile(fileobj=r.raw)
    tar_file = tarfile.open(
        fileobj=io.BytesIO(r.content),
        # mode="r:gz",
        mode="r",
    )
    tar_file.extractall(dir)

    return tar_file


@pyp.utils.memoize
def _map_gene(mapper, symbol_mapping, gene, species):
    symbol = gene
    entrez = brs.mapping.get_entrez_mapping(symbol, species=species)

    if species.lower().replace(" ", "_") not in ["homo_sapiens", "human"]:
        entrez = [i for i in mapper.map_ids([entrez]) if i]

        if not entrez:
            return None, None

        entrez = int(entrez[0])
        symbol = symbol_mapping.loc[entrez]

    return entrez, symbol


def _get_templates(template_dir):
    pyp.utils.makedirs(template_dir)

    url = (
        "https://raw.githubusercontent.com/jdrudolph/photon/"
        "c73d1eb7f5e7cab86031e056350c7b09fb5e1d51/templates/result.html"
    )
    r = requests.get(url)
    r.raise_for_status()

    with open(os.path.join(template_dir, "result.html"), "wb") as f:
        f.write(r.content)


def photon(ds):
    import phos
    import phos.defaults
    import phos.pipeline

    ds = ds.filter(fn=lambda x: len(x["Proteins"]) < 2)
    ds = ds.filter(fn=lambda x: len(x["Scan"]) >= 2)
    ds = ds.filter(mod="Phospho")

    species = list(ds.species)[0]

    from_name = "{}{}".format(
        species.split(" ")[0][0],
        species.split(" ")[1],
    ).lower()

    species = species.replace(" ", "_")

    mapper = EnsemblMapper(
        from_type='entrez',
        to_type='entrez',
        from_organism=from_name,
        to_organism='hsapiens',
    )

    symbol_mapping = brs.cache.get_mapping_data(species="Human")
    symbol_mapping["Symbol"] = symbol_mapping.index
    symbol_mapping = symbol_mapping.set_index("GeneID")["Symbol"]

    def _get_phos_data(psms):
        for _, row in psms.iterrows():
            gene = row["Proteins"].genes[0]

            entrez, symbol = _map_gene(mapper, symbol_mapping, gene, species)

            if not entrez:
                # print(gene, entrez)
                continue

            for mod in row["Modifications"].get_mods("Phospho"):
                yield pd.Series(OrderedDict([
                    ("GeneID", entrez),
                    ("Amino.Acid", mod.letter),
                    ("Position", 1 + mod.abs_pos[0]),
                    ("avg", row["Fold Change"]),
                    ("Symbol", symbol),
                ]))

    LOGGER.info("Generated data frame: {}, {}".format(ds.name, ds.shape))

    df = pd.DataFrame(
        list(_get_phos_data(ds.psms))
    ).dropna()

    LOGGER.info("Generated data frame: {}, {}".format(df.shape))

    df = df.sort_values("avg", ascending=False)

    name = str(uuid.uuid4())

    _parameters = parameters.copy()
    _parameters['anat']['anchor'] = 1950

    dir = os.path.join(pyp.paths.FIGURES_DIR)
    defaults = phos.defaults.make_defaults(dir)
    _get_anat(os.path.abspath("."))

    template_dir = os.path.join(defaults['root'], 'templates')
    _get_templates(template_dir)

    with tempfile.TemporaryDirectory() as work_dir:
        with tempfile.NamedTemporaryFile() as csv_file:
            df.to_csv(csv_file, index=False)

            results = phos.pipeline.run(
                name,
                csv_file.name,
                _parameters,
                template_dir,
                work_dir,
                defaults['db'],
            )

    return results
