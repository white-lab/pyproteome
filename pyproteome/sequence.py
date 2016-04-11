"""
This module provides functionality for manipulating sequences.
"""

# Built-ins
import logging


LOGGER = logging.getLogger("pyproteome.sequence")


class ProteinMatch:
    """
    Contains information mapping a sequence onto a protein.

    Attributes
    ----------
    protein : :class:`Protein<pyproteome.protein.Protein>`
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
    protein_matches : list of :class:`ProteinMatch<pyproteome.sequence.ProteinMatch>`
    modifications : :class:`Modifications<pyproteome.modification.Modifications>`
    alt_hits : list of :class:`Sequence<pyproteome.sequence.Sequence>`
        This attribute indicates which peptides were identified as
        non-ambiguous subsequences of this peptide.
    """

    def __init__(self, pep_seq, protein_matches, modifications=None):
        """
        Parameters
        ----------
        pep_seq : str
        protein_matches : list of :class:`ProteinMatch<pyproteome.sequence.ProteinMatch>`
        modifications : :class:`Modifications<pyproteome.modification.Modifications>`, optional
        """
        self.pep_seq = pep_seq
        self.protein_matches = protein_matches
        self.modifications = modifications
        self.alt_hits = []

    def to_tuple(self):
        return (
            self.pep_seq,
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

        return self.to_tuple() == other.to_tuple()

    def __contains__(self, other):
        if not isinstance(other, Sequence):
            raise TypeError(type(other))

        self_mods = list(self.modifications.skip_labels_iter())
        other_mods = list(other.modifications.skip_labels_iter())

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
            for mod in self.modifications.skip_labels_iter():
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
