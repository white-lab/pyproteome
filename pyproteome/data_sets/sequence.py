"""
This module provides functionality for manipulating sequences.
"""

# Built-ins
import logging

from . import modification, protein

import pyproteome as pyp


LOGGER = logging.getLogger("pyproteome.sequence")


class ProteinMatch:
    """
    Contains information mapping a sequence onto a protein.

    Attributes
    ----------
    protein : :class:`.protein.Protein`
    rel_pos : int
    exact : bool
    """

    def __init__(self, protein, rel_pos, exact):
        self.protein = protein
        self.rel_pos = rel_pos
        self.exact = exact

    def to_tuple(self):
        return (
            self.protein,
            self.rel_pos,
            self.exact,
        )

    def __hash__(self):
        return hash(self.to_tuple())

    def __lt__(self, other):
        return self.protein < other.protein

    def __eq__(self, other):
        if not isinstance(other, ProteinMatch):
            raise TypeError(other)

        return self.to_tuple() == other.to_tuple()


class Sequence:
    """
    Contains information about a sequence and which proteins it matches to.

    Attributes
    ----------
    pep_seq : str
    protein_matches :
    protein_matches : list of :class:`.ProteinMatch`
    modifications : :class:`.modification.Modifications`
    """

    def __init__(
        self,
        pep_seq="",
        protein_matches=None,
        modifications=None,
    ):
        """
        Parameters
        ----------
        pep_seq : str
        protein_matches : list of :class:`.ProteinMatch`
        modifications : :class:`.modification.Modifications`, optional
        """
        if protein_matches is None:
            protein_matches = ()

        self.pep_seq = pep_seq
        self.protein_matches = tuple(sorted(protein_matches))
        self.modifications = modifications

        self._is_labeled = None
        self._is_underlabeled = None

    def to_tuple(self):
        return (
            self.pep_seq.upper(),
            self.modifications,
        )

    def __hash__(self):
        return hash(self.to_tuple())

    def __eq__(self, other):
        # In case of searching just by sequence
        if isinstance(other, str):
            return self._seq_with_modifications() == other

        if not isinstance(other, Sequence):
            raise TypeError(other)

        if self.pep_seq.upper() != other.pep_seq.upper():
            return False

        if tuple(self.protein_matches) != tuple(other.protein_matches):
            return False

        return self.to_tuple() == other.to_tuple()

    def __lt__(self, other):
        return self.pep_seq < other.pep_seq

    def __contains__(self, other):
        if isinstance(other, str):
            return other in self._seq_with_modifications()

        if not isinstance(other, Sequence):
            raise TypeError(type(other))

        self_mods = list(self.modifications.skip_labels())
        other_mods = list(other.modifications.skip_labels())

        return (
            other.pep_seq.upper() in self.pep_seq.upper() and
            len(other.protein_matches) == len(self.protein_matches) and
            all(
                i.protein == j.protein
                for i, j in zip(other.protein_matches, self.protein_matches)
            ) and
            len(other_mods) == len(self_mods) and
            all(
                i.mod_type == j.mod_type and
                i.abs_pos == j.abs_pos and
                i.nterm == j.nterm and
                i.cterm == j.cterm
                for i, j in zip(
                    other_mods,
                    self_mods,
                )
            )
        )

    def _seq_with_modifications(self):
        string = self.pep_seq.upper()

        if self.modifications:
            for mod in self.modifications.skip_labels():
                string = (
                    string[:mod.rel_pos] +
                    string[mod.rel_pos].lower() +
                    string[mod.rel_pos + 1:]
                )

        return string

    def __str__(self):
        if not self.modifications:
            return self.pep_seq
        else:
            return self._seq_with_modifications()

    def __len__(self):
        return len(self.pep_seq)

    @property
    def is_labeled(self):
        """
        Checks whether a sequence is modified on any residue with a
        quantification label.

        Returns
        -------
        is_labeled : bool
        """
        if self._is_labeled is not None:
            return self._is_labeled

        val = any(
            j.mod_type in modification.LABEL_NAMES
            for j in self.modifications.mods
        )

        self._is_labeled = val

        return val

    @property
    def is_underlabeled(self):
        """
        Checks whether a sequence is modified with quantification labels on
        fewer than all expected residues.

        Returns
        -------
        is_underlabeled : bool
        """

        if self._is_underlabeled is not None:
            return self._is_underlabeled

        underlabeled = False

        if self.is_labeled:
            # XXX: Hardcodes label modification locations, not extendable to
            # new quantification tags without changes to this function
            underlabeled = not any(
                j.mod_type in modification.LABEL_NAMES and j.nterm
                for j in self.modifications.mods
            ) or self.pep_seq.count("K") != sum(
                j.mod_type in modification.LABEL_NAMES
                for j in self.modifications.mods
                if j.letter == "K" and not j.nterm
            )

        self._is_underlabeled = underlabeled

        return underlabeled


def extract_sequence(proteins, sequence_string):
    """
    Extract a Sequence object from a list of proteins and sequence string.

    Does not set the Sequence.modifications attribute.

    Parameters
    ----------
    proteins : list of :class:`.protein.Protein`
    sequence_string : str

    Returns
    -------
    seqs : list of :class:`.Sequence`
    """
    prot_matches = []

    # Skip peptides with no protein matches
    if not isinstance(proteins, protein.Proteins):
        proteins = []

    def _get_rel_pos(protein, pep_seq):
        seq = protein.full_sequence

        if not seq:
            return 0, False

        pep_pos = seq.find(pep_seq)
        exact = True

        if pep_pos < 0:
            pep_pos = pyp.utils.fuzzy_find(pep_seq, seq)
            exact = False

        return pep_pos, exact

    for prot in proteins:
        rel_pos, exact = _get_rel_pos(prot, sequence_string.upper())

        prot_matches.append(
            ProteinMatch(
                protein=prot,
                rel_pos=rel_pos,
                exact=exact,
            )
        )

    return Sequence(
        pep_seq=sequence_string,
        protein_matches=prot_matches,
    )
