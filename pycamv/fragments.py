"""
This module provides functionality for calculating the masses of peptide
fragments.
"""

from collections import defaultdict

from . import masses, ms_labels


def _sequence_mass(pep_seq):
    return sum(
        masses.AMINO_ACIDS[letter] +
        sum(
            masses.MODIFICATIONS[letter, mod]
            for mod in mods
        )
        for letter, mods in pep_seq
    )


def _sequence_name(pep_seq):
    return "".join(
        letter.lower() if mods else letter.upper()
        for letter, mods in pep_seq
        if letter not in ["N-term", "C-term"]
    )


def internal_fragment_ions(
    pep_seq,
    any_losses=None, aa_losses=None, mod_losses=None,
):
    """
    Calculate the m/z of all internal fragmenets of a peptide.

    Parameters
    ----------
    pep_seq : list of tuple of (str, list of str)
        A list of peptide letters and each residue's modification(s).
    aa_losses : list of str, optional
        Potential neutral losses for each fragment (i.e. Water, amine, CO).
        List is composed of neutral loss names.
    aa_losses : dict of str, str, optional
    mod_losses : dict of tuple of (str, str), list of str
        Potential neutral losses for modified amino acids (i.e. pY-HPO_3).
        Dictionary should map (letter, modification) to a list of neutral
        loss names.

    Returns
    -------
    dict of str, float
        Dictionary mapping ion names to ion m/z's.
    """
    if any_losses is None:
        any_losses = [
            "H_2O",
            "NH_3",
            "CO",
        ]

    if aa_losses is None:
        aa_losses = {}

    if mod_losses is None:
        mod_losses = {
            ("M", "Oxidation"): ["SOCH_4"],
            ("M", "Dioxidation"): ["SO_2CH_4"],
            ("S", "Phospho"): ["H_3PO_4"],
            ("T", "Phospho"): ["H_3PO_4"],
            ("Y", "Phospho"): ["HPO_3", "HPO_3-H_2O"],
        }

    pep_seq = [
        (letter, mods)
        for letter, mods in pep_seq
        if letter not in ["N-term", "C-term"]
    ]

    frag_masses = {}

    # Calculate the mass of the peptide cut at any two sites between the N-
    # and C-terminus
    for start in range(1, len(pep_seq)):
        for end in range(start + 1, len(pep_seq)):
            # Only add the mass of an N-terminus, cleavage will be between
            # C=O and N-H bond, adding a hydrogen to N-H
            fragment = [("N-term", [])] + pep_seq[start:end]

            mass = _sequence_mass(fragment)
            name = _sequence_name(fragment)

            frag_masses[name] = mass

            # Also calculate any potential neutral losses from this fragment
            losses = _generate_losses(
                fragment, max_depth=1,
                any_losses=any_losses,
                aa_losses=aa_losses,
                mod_losses=mod_losses,
            )

            for loss_name, loss_mass in losses:
                full_loss_name = "{}-{}".format(name, loss_name)
                frag_masses[full_loss_name] = mass - loss_mass

    return frag_masses


def _get_frag_masses(pep_seq):
    return [
        _sequence_mass([pep_seq[index]])
        for index in range(len(pep_seq))
    ]


def _charged_m_zs(name, mass, max_charge):
    for charge in range(1, max_charge + 1):
        yield (
            (
                name.split("-")[0] +
                (
                    "^{{{:+}}}".format(charge)
                    if charge > 1 else
                    "^{+}"
                ) +
                "-".join([""] + name.split("-")[1:])
            ),
            (mass + charge * masses.PROTON) / charge,
        )


def _generate_losses(
    pep_seq,
    max_depth=2,
    any_losses=None, aa_losses=None, mod_losses=None,
):
    def _generate_loss_combos(
        seq,
        losses=None, max_depth=2,
    ):
        if losses is None:
            losses = defaultdict(int)

        if max_depth < 1:
            yield losses
            return

        for loss in any_losses:
            new_losses = losses.copy()
            new_losses[loss] += 1

            yield new_losses

            child_losses = _generate_loss_combos(
                seq, new_losses,
                max_depth=max_depth - 1,
            )

            for new_losses in child_losses:
                yield new_losses

        for aa, a_losses in aa_losses.items():
            for index, (letter, mods) in enumerate(seq):
                if aa != letter:
                    continue

                new_seq = seq[:index] + seq[index + 1:]

                for loss in a_losses:
                    new_losses = losses.copy()
                    new_losses[loss] += 1

                    yield new_losses

                    child_losses = _generate_loss_combos(
                        new_seq, new_losses,
                        max_depth=max_depth - 1,
                    )

                    for new_losses in child_losses:
                        yield new_losses

        for (aa, mod), m_losses in mod_losses.items():
            for index, (letter, mods) in enumerate(seq):
                if aa != letter or mod not in mods:
                    continue

                new_seq = seq[:index] + seq[index + 1:]

                for loss in m_losses:
                    new_losses = losses.copy()
                    new_losses[loss] += 1

                    yield new_losses

                    child_losses = _generate_loss_combos(
                        new_seq, new_losses,
                        max_depth=max_depth - 1,
                    )

                    for new_losses in child_losses:
                        yield new_losses

    losses = _generate_loss_combos(
        pep_seq,
        max_depth=max_depth,
    )
    for loss in losses:
        loss_mass = sum(
            masses.MASSES[loss_name] * count
            for loss_name, count in loss.items()
        )
        loss_name = "-".join(
            "{}{}".format(
                "{} ".format(count) if count > 1 else "",
                loss_name,
            )
            for loss_name, count in loss.items()
        )

        yield loss_name, loss_mass


def _b_y_ions(
    pep_seq, frag_masses, fragment_max_charge,
    any_losses=None, aa_losses=None, mod_losses=None,
):
    def _generate_ions(seq, mass, basename):
        # b/y ions
        charged_mzs = _charged_m_zs(basename, mass, fragment_max_charge)

        for name, mz in charged_mzs:
            yield name, mz

        losses = _generate_losses(
            seq,
            any_losses=any_losses,
            aa_losses=aa_losses,
            mod_losses=mod_losses,
        )

        # b/y ions with losses
        for loss_name, loss_mass in losses:
            charged_mzs = _charged_m_zs(
                "{}-{}".format(basename, loss_name),
                mass - loss_mass,
                fragment_max_charge,
            )

            for name, mz in charged_mzs:
                yield name, mz

    for index in range(2, len(pep_seq) - 1):
        # XXX: iTRAQ / TMT y-adducts?
        # y ion: 1 hydrogen added to NH group, one hydrogen on K/R
        base_ions = {
            "a_{{{}}}".format(index - 1):
                sum(frag_masses[:index]) - masses.MASSES["CO"],
            "b_{{{}}}".format(index - 1):
                sum(frag_masses[:index]),
        }

        for ion_name, ion_mass in base_ions.items():
            for name, mz in _generate_ions(
                pep_seq[:index], ion_mass, ion_name,
            ):
                yield name, mz

    for index in range(1, len(pep_seq) - 1):
        base_ions = {
            "y_{{{}}}".format(len(pep_seq) - index - 1):
                sum(frag_masses[index:]) + masses.PROTON,
        }

        for ion_name, ion_mass in base_ions.items():
            for name, mz in _generate_ions(
                pep_seq[:index], ion_mass, ion_name,
            ):
                yield name, mz


def _label_ions(pep_seq):
    label_mods = [
        mod
        for mod in pep_seq[0][1]
        if mod in ms_labels.LABEL_NAMES
    ]

    for mod in label_mods:
        for name, mz in zip(
            ms_labels.LABEL_NAMES[mod],
            ms_labels.LABEL_MASSES[mod],
        ):
            yield name, mz


def _parent_ions(
    pep_seq, frag_masses, parent_max_charge,
    any_losses=None, aa_losses=None, mod_losses=None,
):
    parent_mass = sum(frag_masses) + masses.PROTON

    for name, mz in _charged_m_zs(
        "MH",
        parent_mass,
        parent_max_charge,
    ):
        yield name, mz

    losses = _generate_losses(
        pep_seq,
        any_losses=any_losses,
        aa_losses=aa_losses,
        mod_losses=mod_losses,
    )

    for loss_name, loss_mass in losses:
        for name, mz in _charged_m_zs(
            "MH",
            parent_mass - loss_mass,
            parent_max_charge,
        ):
            yield "{}-{}".format(name, loss_name), mz


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
    pep_seq, charge,
    parent_max_charge=None, fragment_max_charge=None,
    any_losses=None, aa_losses=None, mod_losses=None,
):
    """
    Calculate the m/z of all ions generated by a peptide.

    Parameters
    ----------
    pep_seq : list of tuple of (str, list of str)
    charge : int
    parent_max_charge : int, optional
    fragment_max_charge : int, optional
    any_losses : list of str, optional
        Potential neutral losses for each fragment (i.e. Water, amine, CO).
        List is composed of neutral loss names.
    aa_losses : dict of str, str, optional
    mod_losses : dict of tuple of (str, str), str, optional
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
        parent_max_charge = charge

    if fragment_max_charge is None:
        # This correct?
        fragment_max_charge = parent_max_charge - 1

    if any_losses is None:
        any_losses = [
            "H_2O",
            "NH_3",
        ]

    if aa_losses is None:
        aa_losses = {}

    if mod_losses is None:
        mod_losses = {
            ("S", "Phospho"): ["H_3PO_4"],
            ("T", "Phospho"): ["H_3PO_4"],
            ("Y", "Phospho"): ["HPO_3", "HPO_3-H_2O"],
            ("M", "Oxidation"): ["SOCH_4"],
            ("M", "Dioxidation"): ["SO_2CH_4"],
        }

    # First calculate the masses of each residue along the backbone
    frag_masses = _get_frag_masses(pep_seq)

    frag_ions = {}

    # Get b/y (and associated a/c/x/z) ions
    frag_ions.update(
        _b_y_ions(
            pep_seq, frag_masses, fragment_max_charge,
            any_losses=any_losses,
            aa_losses=aa_losses,
            mod_losses=mod_losses,
        )
    )

    # Get parent ions (i.e. MH^{+1})
    frag_ions.update(
        _parent_ions(
            pep_seq, frag_masses, parent_max_charge,
            any_losses=any_losses,
            aa_losses=aa_losses,
            mod_losses=mod_losses,
        )
    )

    # Get TMT / iTRAQ labels
    frag_ions.update(
        _label_ions(pep_seq)
    )

    # Get pY peak
    frag_ions.update(
        _py_ions(pep_seq)
    )

    # Get internal fragments
    frag_ions.update(
        internal_fragment_ions(
            pep_seq,
            aa_losses=aa_losses,
            mod_losses=mod_losses,
        )
    )

    return frag_ions
