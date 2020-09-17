'''
This module includes functions for mapping spcies names.
'''


ORGANISM_MAPPING = {
    # 'cat': ,
    # 'chicken': ,
    'Bos taurus': 'cow',
    'Canis familiaris': 'dog',
    'Mustela putorius': 'ferret',
    # 'frog': ,
    'Drosophila melanogaster': 'fruit fly',
    # 'goat': ,
    # 'guinea pig': ,
    # 'hamster': ,
    'Equus caballus': 'horse',
    'Homo sapiens': 'human',
    # 'monkey': ,
    'Mus musculus': 'mouse',
    # 'papillomavirus': ,
    # 'pig': ,
    # 'quail': ,
    # 'rabbit': ,
    'Rattus norvegicus': 'rat',
    # 'sheep': ,
    # 'starfish': ,
    # 'torpedo': ,
    # 'turkey': ,
    # 'water buffalo': ,
}
'''
Mapping between species' specific name and its colloquial name.
(i.e. 'Homo sapiens' > 'human')
'''
INV_ORGANISM_MAPPING = {
    val: key
    for key, val in ORGANISM_MAPPING.items()
}
'''
Mapping between species' colloquial name and its specific name.
'''
