"""
This module provides functionality for calculating the masses of peptide
fragments.
"""

from . import masses


def internal_fragment_masses(pep_seq, aa_losses=None, mod_losses=None):
    """
    Calculate the mass of all internal fragmenets of a peptide.

    Parameters
    ----------
    pep_seq : list of tuple of (str, list of str)
        A list of peptide letters and each residue's modification(s).
    aa_losses : list of str, optional
        Potential neutral losses for each fragment (i.e. Water, amine, CO).
        List is composed of neutral loss names.
    mod_losses : dict of tuple of (str, str), list of str
        Potential neutral losses for modified amino acids (i.e. pY-HPO_3).
        Dictionary should map (letter, modification) to a list of neutral
        loss names.

    Returns
    -------
    dict of str, float
    """
    if aa_losses is None:
        aa_losses = [
            "-H_2O",
            "-NH_3",
            "-CO",
        ]

    if mod_losses is None:
        mod_losses = {
            ("M", "Oxidation"): ["-SOCH_4"],
            ("M", "Dioxidation"): ["-SO_2CH_4"],
            ("S", "Phospho"): ["-H_3PO_4"],
            ("T", "Phospho"): ["-H_3PO_4"],
            ("Y", "Phospho"): ["-HPO_3", "-HPO_3-H_2O"],
        }

    frag_masses = {}

    for start in range(2, len(pep_seq)):
        for end in range(start + 1, len(pep_seq)):
            fragment = pep_seq[start:end]

            # Only add the mass of an N-terminus, cleavage will be between
            # C=O and N-H bond, adding a hydrogen to N-H
            n_term = masses.MASSES[("N-term", None)]

            # XXX: Support multiple modifications to a single residue
            mass = n_term + sum(
                masses.MASSES[(letter, mods[0] if mods else None)]
                for letter, mods in fragment
            )
            name = "".join(
                [
                    letter.lower() if mod else letter
                    for letter, mod in fragment
                ]
            )

            frag_masses[name] = mass

            for loss in aa_losses:
                frag_masses[name + loss] = \
                    mass - masses.MASSES[(loss, None)]

            # M, ST, Y losses
            for (letter, mod), losses in mod_losses.items():
                if not any(
                    letter == l and mod in mods
                    for l, mods in fragment
                ):
                    continue

                for loss in losses:
                    frag_masses[name + loss] = \
                        mass - masses.MASSES[(loss, None)]

    return frag_masses


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
