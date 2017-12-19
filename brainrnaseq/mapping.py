
from . import cache


def get_mapping(gene, species="Mouse"):
    data = cache.get_mapping_data(species=species)

    rows = data[
        data["Symbol"] == gene
    ]

    if rows.shape[0] > 0:
        return rows.iloc[0]["Symbol"]

    rows = data[
        data["Synonyms"].apply(lambda x: gene in x.split("|"))
    ]
    if rows.shape[0] > 0:
        return rows.iloc[0]["Symbol"]

    return None
