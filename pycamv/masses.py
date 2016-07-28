"""
This module provides functionality for calculating molecular masses.
"""

# Built-ins
from __future__ import division

from collections import Iterable

# Core data analysis libraries


def exact_mass(atoms):
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
    >>> exact_mass({"H": 1})        # One hydrogen (Average isotope)
    1.007940721855
    >>> exact_mass({"H": [1]})      # One hydrogen (Exact isotope)
    1.007825
    >>> exact_mass({"H": [1, 1]})   # One hydrogen, one deuterium
    3.021927
    >>> exact_mass({"H": [0, 1]})   # One deuterium (Exact isotope)
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

PROTON = exact_mass({"H": [1]})

C_TERM = exact_mass({"H": [1], "O": [1]})
N_TERM = exact_mass({"H": [1]})

# iTRAQ Masses not exact due to mass differences between isotopes
ITRAQ_4_PLEX = exact_mass({"C": [4, 3], "N": [1, 1], "O": [1], "H": [12]})
ITRAQ_8_PLEX = exact_mass({"C": [7, 7], "N": [3, 1], "O": [3], "H": [24]})
TMT_6_PLEX = exact_mass({"C": [8, 4], "N": [1, 1], "O": [2], "H": [20]})
TMT_10_PLEX = TMT_6_PLEX

LYSINE_HYDROGEN = exact_mass({"H": [1]})

ALANINE = exact_mass({"H": [5], "C": [3], "O": [1], "N": [1]})
ARGININE = exact_mass({"H": [12], "C": [6], "O": [1], "N": [4]})
ASPARAGINE = exact_mass({"H": [6], "C": [4], "O": [2], "N": [2]})
ASPARTATE = exact_mass({"H": [5], "C": [4], "O": [3], "N": [1]})
CYSTEINE = exact_mass({"H": [5], "C": [3], "O": [1], "N": [1], "S": [1]})
GLUTAMATE = exact_mass({"H": [7], "C": [5], "O": [3], "N": [1]})
GLUTAMINE = exact_mass({"H": [8], "C": [5], "O": [2], "N": [2]})
GLYCINE = exact_mass({"H": [3], "C": [2], "O": [1], "N": [1]})
HISTIDINE = exact_mass({"H": [7], "C": [6], "O": [1], "N": [3]})
ISOLEUCINE = exact_mass({"H": [11], "C": [6], "O": [1], "N": [1]})
LEUCINE = exact_mass({"H": [11], "C": [6], "O": [1], "N": [1]})
LYSINE = exact_mass({"H": [12], "C": [6], "O": [1], "N": [2]})
METHIONINE = exact_mass({"H": [9], "C": [5], "O": [1], "N": [1], "S": [1]})
PHENYLALANINE = exact_mass({"H": [9], "C": [9], "O": [1], "N": [1]})
PROLINE = exact_mass({"H": [7], "C": [5], "O": [1], "N": [1]})
SERINE = exact_mass({"H": [5], "C": [3], "O": [2], "N": [1]})
THREONINE = exact_mass({"H": [7], "C": [4], "O": [2], "N": [1]})
TYROSINE = exact_mass({"H": [9], "C": [9], "O": [2], "N": [1]})
TRYPTOPHAN = exact_mass({"H": [10], "C": [11], "O": [1], "N": [2]})
VALINE = exact_mass({"H": [9], "C": [5], "O": [1], "N": [1]})

AMINE = exact_mass({"N": [1], "H": [3]})
WATER = exact_mass({"H": [2], "O": [1]})
OXIDE = exact_mass({"O": [1]})
DIOXIDE = exact_mass({"O": [2]})
CARBON_MONOXIDE = exact_mass({"C": [1], "O": [1]})
CARBON_DIOXIDE = exact_mass({"C": [1], "O": [2]})
PHOSPHORIC_ACID = exact_mass({"H": [3], "P": [1], "O": [4]})
PHOSPHITE = exact_mass({"H": [1], "P": [1], "O": [3]})
SOCH4 = exact_mass({"S": [1], "O": [1], "C": [1], "H": [4]})
SO2CH4 = exact_mass({"S": [1], "O": [2], "C": [1], "H": [4]})
HYDROXYL = exact_mass({"H": [1], "O": [1]})
H2_PHOSPHATE = exact_mass({"H": [2], "O": [4], "P": [1]})
ACETYL = exact_mass({"H": [3], "C": [2], "O": [1]})
CARBAMIDOMETHYL = exact_mass({"H": [2], "C": [2], "O": [1], "N": [1]})

ACETYL_LYSINE = -PROTON + ACETYL
OXY_METHIONINE = OXIDE
DIOXY_METHIONINE = DIOXIDE
CARBAMIDOMETHYL_CYSTEINE = -PROTON + CARBAMIDOMETHYL
PHOSPHO_SERINE = -HYDROXYL + H2_PHOSPHATE
PHOSPHO_THREONINE = -HYDROXYL + H2_PHOSPHATE
PHOSPHO_TYROSINE = -HYDROXYL + H2_PHOSPHATE
PHOSPHO_TYROSINE_IMMONIUM = PHOSPHO_TYROSINE - CARBON_MONOXIDE + PROTON

SILAC_LYSINE_13C6 = exact_mass({"C": [-6, 6]})
SILAC_LYSINE_13C615N2 = exact_mass({"C": [-6, 6], "N": [-2, 2]})
SILAC_ARGININE_13C6 = exact_mass({"C": [-6, 6]})
SILAC_ARGININE_13C615N4 = exact_mass({"C": [-6, 6], "N": [-4, 4]})

# XXX: SILAC Tyrosine?
# XXX: SILAC Leucine?
# XXX: SILAC Acetyl-lysine? (Multiple modifications)

AMINO_ACIDS = {
    # Amino Acids
    "A": ALANINE,
    "R": ARGININE,
    "N": ASPARAGINE,
    "D": ASPARTATE,
    "C": CYSTEINE,
    "E": GLUTAMATE,
    "Q": GLUTAMINE,
    "G": GLYCINE,
    "H": HISTIDINE,
    "I": ISOLEUCINE,
    "L": LEUCINE,
    "K": LYSINE,
    "M": METHIONINE,
    "F": PHENYLALANINE,
    "P": PROLINE,
    "S": SERINE,
    "T": THREONINE,
    "Y": TYROSINE,
    "W": TRYPTOPHAN,
    "V": VALINE,

    # Unknown Amino Acid
    "X": 0,

    # N- and C-terminus
    "N-term": N_TERM,
    "C-term": C_TERM,
}

IMMONIUM_IONS = {
    residue: mass - CARBON_MONOXIDE + PROTON
    for residue, mass in AMINO_ACIDS.items()
}

MODIFICATIONS = {
    # Modified Amino Acids
    ("K", "Acetyl"): ACETYL_LYSINE,
    ("M", "Oxidation"): OXY_METHIONINE,
    ("N", "Dioxidation"): DIOXY_METHIONINE,
    ("C", "Carbamidomethyl"): CARBAMIDOMETHYL_CYSTEINE,
    ("S", "Phospho"): PHOSPHO_SERINE,
    ("T", "Phospho"): PHOSPHO_THREONINE,
    ("Y", "Phospho"): PHOSPHO_TYROSINE,

    # TMT / iTRAQ
    ("N-term", "iTRAQ4plex"): ITRAQ_4_PLEX - N_TERM,
    ("N-term", "iTRAQ8plex"): ITRAQ_8_PLEX - N_TERM,
    ("N-term", "TMT6plex"): TMT_6_PLEX - N_TERM,
    ("N-term", "TMT10plex"): TMT_10_PLEX - N_TERM,
    ("K", "iTRAQ4plex"): ITRAQ_4_PLEX,
    ("K", "iTRAQ8plex"): ITRAQ_8_PLEX,
    ("K", "TMT6plex"): TMT_6_PLEX,
    ("K", "TMT10plex"): TMT_10_PLEX,

    # SILAC
    ("K", "Lysine-13C6 (K-13C6)"): SILAC_LYSINE_13C6,
    ("K", "Lysine-13C615N2 (K-full)"): SILAC_LYSINE_13C615N2,
    ("R", "Arginine-13C6 (R-13C6)"): SILAC_ARGININE_13C6,
    ("R", "Arginine-13C615N2 (R-full)"): SILAC_ARGININE_13C615N4,
}

# Dictionary mapping amino acids and their modifications to weights
MASSES = {
    "Proton": PROTON,

    # Neural Losses
    "H_2O": WATER,
    "NH_3": AMINE,
    "CO": CARBON_MONOXIDE,
    "CO_2": CARBON_DIOXIDE,
    "SOCH_4": SOCH4,                # oxy-M
    "SO2CH4": SO2CH4,               # dioxy-M
    "H_3PO_4": PHOSPHORIC_ACID,     # pS/T
    "HPO_3": PHOSPHITE,
    "HPO_3-H_2O": PHOSPHORIC_ACID,  # pY
}
