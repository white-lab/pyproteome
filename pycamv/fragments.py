"""
This module provides functionality for calculating the masses of peptide
fragments.
"""

import numpy as np

from . import masses, ms_labels


def _get_default_losses(aa_losses, mod_losses):
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

    return aa_losses, mod_losses


def _sequence_mass(pep_seq):
    return sum(
        masses.AMINO_ACIDS[letter] +
        masses.MODIFICATIONS[(letter, mods[0] if mods else None)]
        for letter, mods in pep_seq
    )


def _sequence_name(pep_seq):
    return "".join(
        letter
        for letter, mods in pep_seq
        if letter not in ["N-term", "C-term"]
    )


def internal_fragment_ions(pep_seq, aa_losses=None, mod_losses=None):
    """
    Calculate the m/z of all internal fragmenets of a peptide.

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
        Dictionary mapping ion names to ion m/z's.
    """
    aa_losses, mod_losses = _get_default_losses(aa_losses, mod_losses)

    pep_seq = [
        (letter, mods)
        for letter, mods in pep_seq
        if letter not in ["N-term", "C-term"]
    ]

    frag_masses = {}

    for start in range(2, len(pep_seq)):
        for end in range(start + 1, len(pep_seq)):
            # Only add the mass of an N-terminus, cleavage will be between
            # C=O and N-H bond, adding a hydrogen to N-H
            fragment = [("N-term", [])] + pep_seq[start:end]

            mass = _sequence_mass(fragment)
            name = _sequence_name(fragment)

            frag_masses[name] = mass

            for loss in aa_losses:
                frag_masses[name + loss] = mass - masses.MASSES[loss]

            # M, ST, Y losses
            for (letter, mod), losses in mod_losses.items():
                if not any(
                    letter == l and mod in mods
                    for l, mods in fragment
                ):
                    continue

                for loss in losses:
                    frag_masses[name + loss] = mass - masses.MASSES[loss]

    return frag_masses


def _get_frag_masses(pep_seq):
    return [
        _sequence_mass([pep_seq[index]])
        for index in range(len(pep_seq))
    ]


def _b_y_ions(
    pep_seq, frag_masses,
    fragment_max_charge,
    aa_losses, mod_losses
):
    proton = masses.PROTON

    ions = {}

    def _generate_ions(mass, basename):
        ret = {}

        ret[basename] = mass

        # iTRAQ / TMT y-adducts?

        for loss in aa_losses + [""]:
            loss_mass = masses.MASSES[loss]

            for charge in range(2, fragment_max_charge):
                ret[
                    name + "^\{{:+}\}".format(charge)
                ] = (mass - loss_mass + charge * proton) / charge

        return ret

    for index in range(2, len(pep_seq) - 1):
        # b-ions first
        mass = np.cumsum(frag_masses[:index])
        name = "b_\{{}\}".format(index - 1)
        ions.update(_generate_ions(mass, name))

        # y-ions second
        mass = np.cumsum(frag_masses[index:])
        name = "y_\{{}\}".format(index - 1)
        ions.update(_generate_ions(mass, name))

    return ions


def _label_ions(pep_seq):
    ions = {}

    label_mods = [
        mod
        for mod in pep_seq[0][1]
        if mod in ms_labels.LABEL_NAMES
    ]

    for mod in label_mods:
        for name, m_z in zip(
            ms_labels.LABEL_NAMES[mod],
            ms_labels.LABEL_MASSES[mod],
        ):
            ions[name] = m_z

    return ions


def _parent_ions(frag_masses, parent_max_charge):
    proton = masses.PROTON
    ions = {}
    parent = sum(frag_masses)

    for charge in range(1, parent_max_charge):
        name = "MH^\{{:+}\}".format(charge)
        ions[name] = (parent + charge * proton) / charge

    return ions


def _py_ions(pep_seq):
    ions = {}

    if any(
        letter == "Y" and "Phospho" in mods
        for letter, mods in pep_seq
    ):
        ions["pY"] = masses.IMMONIUM_IONS["Y"] + \
            masses.MODIFICATIONS["Y", "Phospho"]

    return ions


def fragment_ions(
    pep_seq,
    parent_max_charge=None, fragment_max_charge=None,
    aa_losses=None, mod_losses=None,
):
    """
    Calculate the m/z of all ions generated by a peptide.

    Parameters
    ----------
    pep_seq : str
    parent_max_charge : int, optional
    fragment_max_charge : int, optional
    aa_losses : list of str, optional
        Potential neutral losses for each fragment (i.e. Water, amine, CO).
        List is composed of neutral loss names.
    mod_losses : dict of tuple of (str, str), list of str
        Potential neutral losses for modified amino acids (i.e. pY-HPO_3).
        Dictionary should map (letter, modification) to a list of neutral
        loss names.

    Returns
    -------
    dict of int, dict of str, float
        Dictionary mapping fragment position to a dictionary mapping ion names
        to ion m/z's.
    """
    assert "N-term" == pep_seq[0][0]
    assert "C-term" == pep_seq[-1][0]

    if parent_max_charge is None:
        parent_max_charge = 5

    if fragment_max_charge is None:
        fragment_max_charge = 4

    fragment_charge_range = range(1, fragment_max_charge)
    aa_losses, mod_losses = _get_default_losses(aa_losses, mod_losses)

    # First calculate the masses of each residue along the backbone
    frag_masses = _get_frag_masses(pep_seq)

    other_ions = {}

    # Get b/y (and associated a/c/x/z) ions
    frag_ions = _b_y_ions(
        pep_seq, frag_masses, fragment_charge_range,
        aa_losses, mod_losses,
    )

    # Get parent ions (i.e. MH^{+1})
    other_ions.update(
        _parent_ions(frag_masses, parent_max_charge)
    )

    # Get TMT / iTRAQ labels
    other_ions.update(
        _label_ions(pep_seq)
    )

    # Get pY peak
    other_ions.update(
        _py_ions(pep_seq)
    )

    # Get internal fragments
    other_ions.update(
        internal_fragment_ions(
            pep_seq,
            aa_losses=aa_losses,
            mod_losses=mod_losses,
        )
    )

    return frag_ions, other_ions
