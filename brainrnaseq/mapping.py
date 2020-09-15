
from . import cache, utils


def _find_uniprot_gene(gene, data):
    try:
        rows = data.loc[gene]
    except KeyError:
        pass
    else:
        if len(rows.shape) > 1:
            rows = rows.iloc[0]

        return rows

    rows = data[
        data['Synonyms'].apply(lambda x: gene in x.split('|'))
    ]

    if rows.shape[0] > 0:
        return rows.iloc[0]

    return None


@utils.memoize
def get_symbol_mapping(gene, species='Mouse'):
    '''
    Returns
    -------
    pandas.Series
    '''
    data = cache.get_mapping_data(species=species)

    row = _find_uniprot_gene(gene, data)
    return row.name if row is not None else None


@utils.memoize
def get_entrez_mapping(gene, species='Mouse'):
    '''
    Returns
    -------
    pandas.Series
    '''
    data = cache.get_mapping_data(species=species)

    row = _find_uniprot_gene(gene, data)
    return str(row['GeneID']) if row is not None else None
