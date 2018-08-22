
import logging
import requests

import pandas as pd

import pyproteome as pyp


MSIGDB_URL = (
    "https://raw.githubusercontent.com/white-lab/pyproteome-data"
    "/master/msigdb/"
)
MSIGDB_FILES = (
    "h.all.v6.1.entrez.gmt",
    # "c1.all.v6.1.entrez.gmt",
    # "c2.all.v6.1.entrez.gmt",
    # "c2.cgp.v6.1.entrez.gmt",
    # "c2.cp.biocarta.v6.1.entrez.gmt",
    # "c2.cp.kegg.v6.1.entrez.gmt",
    # "c2.cp.reactome.v6.1.entrez.gmt",
    # "c2.cp.v6.1.entrez.gmt",
    # "c3.all.v6.1.entrez.gmt",
    # "c3.mir.v6.1.entrez.gmt",
    # "c3.tft.v6.1.entrez.gmt",
    # "c4.all.v6.1.entrez.gmt",
    # "c4.cgn.v6.1.entrez.gmt",
    # "c4.cm.v6.1.entrez.gmt",
    # "c5.all.v6.1.entrez.gmt",
    # "c5.bp.v6.1.entrez.gmt",
    # "c5.cc.v6.1.entrez.gmt",
    # "c5.mf.v6.1.entrez.gmt",
    # "c6.all.v6.1.entrez.gmt",
    # "c7.all.v6.1.entrez.gmt",
)

LOGGER = logging.getLogger("pyproteome.msigdb")

try:
    from genemap.mappers import EnsemblMapper
except ImportError:
    pass


@pyp.utils.memoize
def get_msigdb_pathways(species, remap=None):
    """
    Download gene sets from MSigDB. Currently downloads v6.1 of the gene
    signature repositories.

    Parameters
    ----------
    species : str

    Returns
    -------
    df : :class:`pandas.DataFrame`, optional
    """
    LOGGER.info("Fetching MSigDB pathways")

    def _get_requests():
        for file in MSIGDB_FILES:
            url = MSIGDB_URL + file

            LOGGER.info("Fetching {}".format(url))

            response = requests.get(url, stream=True)
            response.raise_for_status()

            yield response

    def _get_data(line):
        line = line.decode("utf-8")
        name, _, genes = line.split("\t", 2)
        # name, _, _, spec = name.split("%")
        # assert species == spec
        return name, set(i for i in genes.split("\t"))

    pathways_df = pd.DataFrame(
        data=[
            _get_data(line)
            for response in _get_requests()
            for line in response.iter_lines()
        ],
        columns=["name", "set"],
    )

    if remap and species not in ["Homo sapiens"]:
        to_name = "{}{}".format(
            species.split(" ")[0][0],
            species.split(" ")[1],
        ).lower()

        LOGGER.info("Remapping MSigDB to {} ({})".format(species, to_name))

        mapper = EnsemblMapper(
            from_type='entrez',
            to_type='entrez',
            from_organism='hsapiens',
            to_organism=to_name,
        )
        pathways_df["set"] = pathways_df["set"].apply(
            lambda row: set(mapper.map_ids(row))
        )

    return pathways_df
