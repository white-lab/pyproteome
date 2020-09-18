
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.go_search import GoSearch
from goatools.anno.genetogo_reader import Gene2GoReader
from collections import defaultdict

import importlib
import re

TAXA = {
    'Mus musculus': 10090,
    'Homo sapiens': 9606,
}

def get_go_ids(go_ids, species='Homo sapiens'):
    '''
    Fetch all gene symbols associated with a list of gene ontology term IDs.

    Parameters
    ----------
    go_ids : str or list of str
    species : str, optional

    Returns
    -------
    list of str
    '''
    assert species in TAXA

    if isinstance(go_ids, str):
        go_ids = [go_ids]

    obo_fname = download_go_basic_obo('db/go/go-basic.obo')
    gene2go = download_ncbi_associations('db/go/gene2go')

    taxid = TAXA[species]

    fin_symbols = 'genes_NCBI_{TAXID}_All.py'.format(TAXID=taxid)

    module_name = ''.join(['goatools.test_data.', fin_symbols[:-3]])
    module = importlib.import_module(module_name)
    GeneID2nt = module.GENEID2NT

    go2geneids = Gene2GoReader(
        'db/go/gene2go',
        taxids=[taxid],
    )

    go2items = defaultdict(list)
    for i in go2geneids.taxid2asscs[taxid]:
        go2items[i.GO_ID].append(i.DB_ID)

    srchhelp = GoSearch('db/go/go-basic.obo', go2items=go2items)

    with open('go.log', 'w') as log:
        # Add children GOs
        gos_all = srchhelp.add_children_gos(go_ids)
        
        # Get Entrez GeneIDs for cell cycle GOs
        gene_ids = set()

        for go_items in [
            go_ids,
            gos_all,
        ]:
            gene_ids.update(srchhelp.get_items(go_items))

    genes = []

    for geneid in gene_ids:
        nt = GeneID2nt.get(geneid, None)

        if nt is not None:
            genes.append(nt.Symbol)

    return genes
