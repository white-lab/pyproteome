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
ACETYL_LYSINE = 0
OXY_METHIONINE = 0
DIOXY_METHIONINE = 0
CARBAMIDOMETHYL_CYSTEINE = 0
PHOSPHO_SERINE = 0
PHOSPHO_THREONINE = 0
PHOSPHO_TYROSINE = 0

AMINE = _exact_mass({"N": [1], "H": [3]})
WATER = _exact_mass({"H": [2], "O": [1]})
PHOSPHORIC_ACID = _exact_mass({"H": [3], "P": [1], "O": [4]})
PHOSPHITE = _exact_mass({"H": [1], "P": [1], "O": [3]})
CARBON_DIOXIDE = _exact_mass({"C": [1], "O": [2]})
SOCH4 = _exact_mass({"S": [1], "O": [1], "C": [1], "H": [4]})
SO2CH4 = _exact_mass({"S": [1], "O": [2], "C": [1], "H": [4]})
