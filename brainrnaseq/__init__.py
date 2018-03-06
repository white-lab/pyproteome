
from . import cache, mapping, enrichments


CELL_TYPE_COLS = {
    "Homo sapiens": {
        "Astrocyte": [
            '8yo',
            '13yo', '16yo', '21yo.1', '22yo.1', '35yo', '47yo', '51yo', '53yo',
            '60yo', '63yo - 1', '63yo - 2',
        ],
        "Neuron": [
            "25yo",
        ],
        "OPC": [
            '22yoGC', '63yoGC - 1',
            '63yo GC - 2', '47yoO4', '63yoO4',
        ],
        "New Oligodendrocytes": [
            '22yoGC', '63yoGC - 1',
            '63yo GC - 2', '47yoO4', '63yoO4',
        ],
        "Myelinating Oligodendrocytes": [
            '22yoGC', '63yoGC - 1',
            '63yo GC - 2', '47yoO4', '63yoO4',
        ],
        "Microglia": [
            '45yo', '51yo.1', '63yo',
        ],
        "Endothelia": [
            "13yo.1",
            "47yo.1",
        ],
    },
    "Mus musculus": {
        "Astrocyte": [
            # "FACS - p69",
            # "FACS p70",
            "1 month",
            "4 months",
            "7 months",
            "9 months",
        ],
        "Neuron": [
            "Neuron 3",
            "Neuron 4",
        ],
        "OPC": [
            "Oligodendrocyte precursor cell 3",
            "Oligodendrocyte precursor cell 4",
        ],
        "New Oligodendrocytes": [
            "Newly formed oligodendrocyte 3",
            "Newly formed oligodendrocyte 4",
        ],
        "Myelinating Oligodendrocytes": [
            "Myelinating oligodendrocyte 4",
            "Myelinating oligodenrocyte 5",
        ],
        "Microglia": [
            "Microglia 1",
            "Microglia 2",
        ],
        "Endothelia": [
            "Endo 1",
            "Endo 2",
        ],
    },
}

CELL_TYPES = [
    "Astrocyte",
    "Neuron",
    "OPC",
    "New Oligodendrocytes",
    "Myelinating Oligodendrocytes",
    "Microglia",
    "Endothelia",
]
DEFAULT_CELL_TYPES = CELL_TYPES[:2] + CELL_TYPES[4:]

__all__ = [
    "cache",
    "mapping",
    "enrichments",

    "CELL_TYPE_COLS",
    "CELL_TYPES",
    "DEFAULT_CELL_TYPES",
]
