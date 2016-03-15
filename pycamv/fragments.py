"""
This module provides functionality for calculating the masses of peptide
fragments.
"""


def internal_fragment_masses(pep_seq):
    """
    Calculate the mass of all internal fragmenets of a peptide.

    Parameters
    ----------
    pep_seq : str

    Returns
    -------
    dict of str, float
    """
    pass


def fragment_masses(pep_seq, charge_range=None):
    """
    Calculate the mass of all b/y ions in a peptide.

    Parameters
    ----------
    pep_seq : str
    charge_range : tuple of (int, int), optional

    Returns
    -------
    dict of str, float
    """
    pass
