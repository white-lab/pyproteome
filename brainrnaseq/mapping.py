
from . import cache


def get_mapping(gene, species="Mouse"):
    data = cache.get_mapping_data(species=species)

    try:
        row = data[
            data["Symbol"] == gene
        ].iloc[0]
    except IndexError:
        row = data[
            data["Synonyms"].apply(lambda x: gene in x.split("|"))
        ].iloc[0]

    return row["Symbol"]
