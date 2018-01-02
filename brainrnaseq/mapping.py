
from . import cache


def get_mapping(gene, species="Mouse", force=False):
    data = cache.get_mapping_data(species=species, force=force)

    try:
        data.loc[gene]
        return gene
    except KeyError:
        pass

    rows = data[
        data["Synonyms"].apply(lambda x: gene in x.split("|"))
    ]
    if rows.shape[0] > 0:
        return rows.index[0]

    return None
