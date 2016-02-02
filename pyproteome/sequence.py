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
    protein : pyproteome.Protein
    rel_pos : int
    exact : bool
    """

    def __init__(self, protein, rel_pos, exact):
        self.protein = protein
        self.rel_pos = rel_pos
        self.exact = exact

    def __hash__(self):
        return hash(
            (
                self.protein,
                self.rel_pos,
                self.exact,
            )
        )

    def __eq__(self, other):
        if not isinstance(other, ProteinMatch):
            raise TypeError()

        return (
            self.protein,
            self.rel_pos,
            self.exact,
        ) == (
            other.protein,
            other.rel_pos,
            other.exact,
        )


class Sequence:
    """
    Contains information about a sequence and which proteins it matches to.

    Attributes
    ----------
    pep_seq : str
    protein_matches : list of pyproteome.ProteinMatch
    modifications : pyproteome.Modifications
    """

    def __init__(self, pep_seq, protein_matches, modifications=None):
        """
        Parameters
        ----------
        pep_seq : str
        protein_matches : list of pyproteome.ProteinMatch
        modifications : pyproteome.Modifications, optional
        """
        self.pep_seq = pep_seq
        self.protein_matches = protein_matches
        self.modifications = modifications

    def __hash__(self):
        return hash(
            (
                self.pep_seq,
                self.modifications,
            )
        )

    def __eq__(self, other):
        # In case of searching just by sequence
        if isinstance(other, str):
            return self._seq_with_modifications() == other

        if not isinstance(other, Sequence):
            raise TypeError()

        return (
            self.pep_seq,
            self.modifications,
        ) == (
            other.pep_seq,
            other.modifications,
        )

    def _seq_with_modifications(self):
        string = self.pep_seq.upper()

        for mod in self.modifications.skip_labels_iter():
            string = string[:mod.rel_pos] + string[mod.rel_pos].lower() + \
                string[mod.rel_pos + 1:]

        return string

    def __str__(self):
        if not self.modifications:
            return self.pep_seq
        else:
            return self._seq_with_modifications()
