"""
This module provides functionality for calculating molecular masses.
"""

# Built-ins
from __future__ import division

from collections import Iterable

# Core data analysis libraries


def _exact_mass(atoms):
    """
    Parameters
    ----------
    atoms : dict of str, int or iterable
        Atoms in molecule and their count. If the dictionary's values are
        integers, use average isotope weights. If values are iterables, use
        exact count of each isotope, in order of increasing weight.

    Returns
    -------
    float
        Molecular weight of molecule.

    Examples
    --------
    >>> _exact_mass({"H": 1})        # One hydrogen (Average isotope)
    1.007940721855
    >>> _exact_mass({"H": [1]})      # One hydrogen (Exact isotope)
    1.007825
    >>> _exact_mass({"H": [1, 1]})   # One hydrogen, one deuterium
    3.021927
    >>> _exact_mass({"H": [0, 1]})   # One deuterium (Exact isotope)
    2.014102
    """
    # Atom name mapping to molecular weight and frequency of that isotope
    molecular_weights = {
        "H": [
            (1.007825, 99.9885),
            (2.014102,  0.0115),
        ],
        "C": [
            (12.000000, 98.93),
            (13.003355,  1.07)
        ],
        "N": [
            (14.003074, 99.632),
            (15.000109,  0.368),
        ],
        "O": [
            (15.994915, 99.757),
            (16.999131,  0.038),
            (17.999159,  0.205),
        ],
        "S": [
            (31.972072, 94.93),
            (32.971459,  0.76),
            (33.967868,  4.29),
            (35.967079,  0.02),
        ],
        "P": [
            (30.973763, 100),
        ],
    }
    return sum(
        sum(
            weight[0] * count
            for weight, count in zip(molecular_weights[atom], counts)
        )
        if isinstance(counts, Iterable) else
        sum(
            weight * counts * frequency / 100
            for weight, frequency in molecular_weights[atom]
        )
        for atom, counts in atoms.items()
    )


C_TERM = _exact_mass({"H": [1], "O": [1]})
N_TERM = _exact_mass({"H": [1]})

# iTRAQ Masses not exact due to mass differences between isotopes
ITRAQ_4_PLEX = _exact_mass({"C": [4, 3], "N": [1, 1], "O": [1], "H": [12]})
ITRAQ_8_PLEX = _exact_mass({"C": [7, 7], "N": [3, 1], "O": [3], "H": [24]})
TMT_6_PLEX = _exact_mass({"C": [8, 4], "N": [1, 1], "O": [2], "H": [20]})
TMT_10_PLEX = TMT_6_PLEX

LYSINE_HYDROGEN = _exact_mass({"H": [1]})

ALANINE = _exact_mass({"H": [5], "C": [3], "O": [1], "N": [1]})
ARGININE = _exact_mass({"H": [12], "C": [6], "O": [1], "N": [4]})
ASPARAGINE = _exact_mass({"H": [6], "C": [4], "O": [2], "N": [2]})
ASPARTATE = _exact_mass({"H": [5], "C": [4], "O": [3], "N": [1]})
CYSTEINE = _exact_mass({"H": [5], "C": [3], "O": [1], "N": [1], "S": [1]})
GLUTAMATE = _exact_mass({"H": [7], "C": [5], "O": [3], "N": [1]})
GLUTAMINE = _exact_mass({"H": [8], "C": [5], "O": [2], "N": [2]})
GLYCINE = _exact_mass({"H": [3], "C": [2], "O": [1], "N": [1]})
HISTIDINE = _exact_mass({"H": [7], "C": [6], "O": [1], "N": [3]})
ISOLEUCINE = _exact_mass({"H": [11], "C": [6], "O": [1], "N": [1]})
LEUCINE = _exact_mass({"H": [11], "C": [6], "O": [1], "N": [1]})
LYSINE = _exact_mass({"H": [12], "C": [6], "O": [1], "N": [2]})
METHIONINE = _exact_mass({"H": [9], "C": [5], "O": [1], "N": [1], "S": [1]})
PHENYLALANINE = _exact_mass({"H": [9], "C": [9], "O": [1], "N": [1]})
PROLINE = _exact_mass({"H": [7], "C": [5], "O": [1], "N": [1]})
SERINE = _exact_mass({"H": [5], "C": [3], "O": [2], "N": [1]})
THREONINE = _exact_mass({"H": [7], "C": [4], "O": [2], "N": [1]})
TYROSINE = _exact_mass({"H": [9], "C": [9], "O": [2], "N": [1]})
TRYPTOPHAN = _exact_mass({"H": [10], "C": [11], "O": [1], "N": [2]})
VALINE = _exact_mass({"H": [9], "C": [5], "O": [1], "N": [1]})

R = 0
ACETYL_LYSINE = \
    _exact_mass({"H": [14], "C": [8], "O": [2], "N": [2]})
OXY_METHIONINE = \
    _exact_mass({"H": [9], "C": [5], "O": [2], "N": [1], "S": [1]})
DIOXY_METHIONINE = \
    _exact_mass({"H": [9], "C": [5], "O": [3], "N": [1], "S": [1]})
CARBAMIDOMETHYL_CYSTEINE = \
    _exact_mass({"H": [8], "C": [5], "O": [2], "N": [2], "S": [1]})
PHOSPHO_SERINE = \
    _exact_mass({"H": [6], "C": [3], "O": [5], "N": [1], "P": [1]})
PHOSPHO_THREONINE = \
    _exact_mass({"H": [8], "C": [4], "O": [5], "N": [1], "P": [1]})
PHOSPHO_TYROSINE = \
    _exact_mass({"H": [10], "C": [9], "O": [5], "N": [1], "P": [1]})

SILAC_LYSINE_13C6 = \
    _exact_mass({"H": [12], "C": [0, 6], "O": [1], "N": [2]})
SILAC_LYSINE_13C615N2 = \
    _exact_mass({"H": [12], "C": [0, 6], "O": [1], "N": [0, 2]})
SILAC_ARGININE_13C6 = \
    _exact_mass({"H": [12], "C": [0, 6], "O": [1], "N": [4]})
SILAC_ARGININE_13C615N4 = \
    _exact_mass({"H": [12], "C": [0, 6], "O": [1], "N": [0, 4]})

# XXX: SILAC Tyrosine?
# XXX: SILAC Leucine?
# XXX: SILAC Acetyl-lysine? (Multiple modifications)

AMINE = _exact_mass({"N": [1], "H": [3]})
WATER = _exact_mass({"H": [2], "O": [1]})
PHOSPHORIC_ACID = _exact_mass({"H": [3], "P": [1], "O": [4]})
PHOSPHITE = _exact_mass({"H": [1], "P": [1], "O": [3]})
CARBON_DIOXIDE = _exact_mass({"C": [1], "O": [2]})
SOCH4 = _exact_mass({"S": [1], "O": [1], "C": [1], "H": [4]})
SO2CH4 = _exact_mass({"S": [1], "O": [2], "C": [1], "H": [4]})


# Dictionary mapping amino acids and their modifications to weights
MASSES = dict(
    [
        (("A", None), ALANINE),
        (("R", None), ARGININE),
        (("N", None), ASPARAGINE),
        (("D", None), ASPARTATE),
        (("C", None), CYSTEINE),
        (("E", None), GLUTAMATE),
        (("Q", None), GLUTAMINE),
        (("G", None), GLYCINE),
        (("H", None), HISTIDINE),
        (("I", None), ISOLEUCINE),
        (("L", None), LEUCINE),
        (("K", None), LYSINE),
        (("M", None), METHIONINE),
        (("F", None), PHENYLALANINE),
        (("P", None), PROLINE),
        (("S", None), SERINE),
        (("T", None), THREONINE),
        (("Y", None), TYROSINE),
        (("W", None), TRYPTOPHAN),
        (("V", None), VALINE),
        (("K", "Acetyl"), ACETYL_LYSINE),
        (("M", "Oxidation"), OXY_METHIONINE),
        (("N", "Dioxidation"), DIOXY_METHIONINE),
        (("C", "Carbamidomethyl"), CARBAMIDOMETHYL_CYSTEINE),
        (("S", "Phospho"), PHOSPHO_SERINE),
        (("T", "Phospho"), PHOSPHO_THREONINE),
        (("Y", "Phospho"), PHOSPHO_TYROSINE),
        (("N-term", None), N_TERM),
        (("C-term", None), C_TERM),
        (("N-term", "iTRAQ4plex"), N_TERM + ITRAQ_4_PLEX),
        (("N-term", "iTRAQ8plex"), N_TERM + ITRAQ_8_PLEX),
        (("N-term", "TMT6plex"), N_TERM + TMT_6_PLEX),
        (("N-term", "TMT10plex"), N_TERM + TMT_10_PLEX),
        (("K", "iTRAQ4plex"), LYSINE + ITRAQ_4_PLEX - LYSINE_HYDROGEN),
        (("K", "iTRAQ8plex"), LYSINE + ITRAQ_8_PLEX - LYSINE_HYDROGEN),
        (("K", "TMT6plex"), LYSINE + TMT_6_PLEX - LYSINE_HYDROGEN),
        (("K", "TMT10plex"), LYSINE + TMT_10_PLEX - LYSINE_HYDROGEN),
        (("K", "Lysine-13C6 (K-13C6)"), SILAC_LYSINE_13C6),
        (("K", "Lysine-13C615N2 (K-full)"), SILAC_LYSINE_13C615N2),
        (("R", "Arginine-13C6 (R-13C6)"), SILAC_ARGININE_13C6),
        (("R", "Arginine-13C615N2 (R-full)"), SILAC_ARGININE_13C615N4),
    ]
)
