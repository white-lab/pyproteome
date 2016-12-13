"""
This module provides functionality for post-translational modifications.

Wraps modifications in a structured class and allows filtering of
modifications by amino acid and modification type.
"""

LABEL_NAMES = ["TMT", "ITRAQ"]


class Modifications:
    """
    A list of modifications.

    Wraps the Modification objects and provides several utility functions.

    Attributes
    ----------
    mods : list of :class:`Modification<pyproteome.modification.Modification>`
    """

    def __init__(self, mods):
        """
        Initialize from a list of modifications.

        Parameters
        ----------
        mods : list of :class:`Modification<pyproteome.modification.Modification>`
        """
        self.mods = mods

    def __iter__(self):
        return iter(self.mods)

    def __hash__(self):
        return hash(
            tuple(self.mods),
        )

    def __len__(self):
        return len(self.mods)

    def skip_labels_iter(self):
        """
        Return an iterable, skipping over any modifications for peptide labels.

        Returns
        -------
        generator of :class:`Modification<pyproteome.modification.Modification>`
        """
        return (
            mod
            for mod in self.mods
            if not any(label in mod.mod_type for label in LABEL_NAMES)
        )

    def get_mods(self, letter_mod_types):
        """
        Filter the list of modifications.

        Only keeps modifications with a given letter, mod_type, or both.

        Parameters
        ----------
        letter_mod_types : list of tuple of str, str

        Returns
        -------
        generator of Modification
        """
        any_letter, any_mod, letter_mod = \
            _extract_letter_mods(letter_mod_types)
        return Modifications(
            [
                mod
                for mod in self.mods
                if allowed_mod_type(
                    mod,
                    any_letter=any_letter,
                    any_mod=any_mod,
                    letter_mod=letter_mod,
                )
            ]
        )

    def __eq__(self, other):
        if not isinstance(other, Modifications):
            raise TypeError()

        # XXX: Incorrect if modifications are sorted differently
        return len(self.mods) == len(other.mods) and all(
            i == j
            for i, j in zip(self.mods, other.mods)
        )

    def __str__(self, absolute=True, skip_labels=True):
        if len(self.mods) == 0:
            return ""

        if skip_labels:
            lst = list(self.skip_labels_iter())
        else:
            lst = list(iter(self))

        if not lst:
            return ""

        return " / ".join(
            ", ".join(
                "{}{}{}{}".format(
                    mod.display_mod_type(),
                    mod.letter.upper(),
                    1 + (mod.abs_pos[i] if absolute else mod.rel_pos),
                    "" if mod.exact[i] else "*"
                )
                for mod in lst
            )
            for i in range(len(lst[0].exact))
        )


class Modification:
    """
    Contains information for a single peptide modification.

    Attributes
    ----------
    rel_pos : int
    mod_type : str
    nterm : bool
    cterm : bool
    letter : str
    abs_pos : int
    """

    def __init__(self, rel_pos, mod_type, sequence, nterm=False, cterm=False):
        self.rel_pos = rel_pos
        self.mod_type = mod_type
        self.nterm = nterm
        self.cterm = cterm
        self.letter = sequence.pep_seq[rel_pos]
        self.abs_pos = [
            rel_pos + match.rel_pos
            for match in sequence.protein_matches
        ]
        self.exact = [
            match.exact
            for match in sequence.protein_matches
        ]

    def display_mod_type(self):
        """
        Return the mod_type in an abbreviated form (i.e. "p" for "Phospho")

        Returns
        -------
        str
        """
        if self.mod_type in ["Phospho"]:
            return "p"
        if self.mod_type in ["Carbamidomethyl"]:
            return "cm"

        return self.mod_type

    def to_tuple(self):
        return (
            self.rel_pos,
            self.mod_type,
            self.nterm,
            self.cterm,
            self.letter,
            tuple(self.abs_pos),
            tuple(self.exact),
        )

    def __hash__(self):
        return hash(self.to_tuple())

    def __eq__(self, other):
        if not isinstance(other, Modification):
            raise TypeError()

        return self.to_tuple() == other.to_tuple()


def allowed_mod_type(mod, any_letter=None, any_mod=None, letter_mod=None):
    """
    Check if a modification is of a type.

    Filters by letter, mod_type, or both.

    Parameters
    ----------
    mod : :class:`Modification<pyproteome.modification.Modification>`
    any_letter : set of str
    any_mod : set of str
    letter_mod : set of tuple of str, str
    """
    return (
        (
            any_letter is None or
            mod.letter.upper() in any_mod
        ) or (
            any_mod is None or
            mod.mod_type in any_letter
        ) or (
            letter_mod is None or
            (mod.letter.upper(), mod.mod_type) in letter_mod
        )
    )


def _extract_letter_mods(letter_mod_types=None):
    if letter_mod_types is None:
        return None, None, None

    any_letter = set()
    any_mod = set()
    letter_mod = set()

    for letter, mod_type in letter_mod_types:
        if letter is None and mod_type is None:
            raise Exception("Need at least one letter or mod type not None")
        elif letter is None and mod_type is not None:
            any_letter.add(mod_type)
        elif letter is not None and mod_type is None:
            any_mod.add(letter.upper())
        else:
            letter_mod.add((letter.upper(), mod_type))

    return any_letter, any_mod, letter_mod


def filter_mod_types(psms, letter_mod_types=None):
    """
    Filter a list of psms sequences with a given modification type.

    Parameters
    ----------
    psms : :class:`pandas.DataFrame`
    letter_mod_types : list of tuple of str or None, str or None, optional

    Returns
    -------
    :class:`pandas.DataFrame`
    """
    return psms[
        psms["Modifications"].apply(
            lambda x: bool(x.get_mods(letter_mod_types).mods)
        )
    ]
