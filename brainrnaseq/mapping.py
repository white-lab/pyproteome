
from . import cache


def _find_uniprot_gene(gene, data):
    try:
        row = data.loc[gene]
        return row
    except KeyError:
        pass

    rows = data[
        data["Synonyms"].apply(lambda x: gene in x.split("|"))
    ]
    if rows.shape[0] > 0:
        return rows.iloc[0]

    return None


def get_symbol_mapping(gene, species="Mouse", force=False):
    """
    Returns
    -------
    pandas.Series
    """
    data = cache.get_mapping_data(species=species, force=force)

    row = _find_uniprot_gene(gene, data)
    return row.index if row is not None else None


def get_entrez_mapping(gene, species="Mouse", force=False):
    """
    Returns
    -------
    pandas.Series
    """
    data = cache.get_mapping_data(species=species, force=force)

    row = _find_uniprot_gene(gene, data)
    return row["GeneID"] if row is not None else None
