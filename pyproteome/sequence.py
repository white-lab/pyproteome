"""
This module provides functionality for manipulating sequences.

Functionality includes n-mer generation.
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


def generate_n_mers(
    sequences, n=15,
    all_matches=True,
    fill_left="A",
    fill_right="A",
    letter_mod_types=None,
):
    """
    Generate n-mers around all sites of modification in sequences.

    Parameters
    ----------
    sequences : list of pyproteome.Sequence
    n : int, optional
    all_matches : bool, optional
    fill_left : str, optional
    fill_right : str, optional
    """
    # Check n is odd
    assert n % 2 == 1

    def _n_mer_from_sequence(full_seq, abs_pos):
        return (
            fill_left * (n // 2 - abs_pos) +
            full_seq[max([abs_pos - n // 2, 0]):abs_pos] +
            full_seq[abs_pos].lower() +
            full_seq[abs_pos + 1:abs_pos + n // 2 + 1] +
            fill_right * (abs_pos - len(full_seq) + n // 2 + 1)
        )

    return (
        _n_mer_from_sequence(
            match.protein.full_sequence,
            mod.rel_pos + match.rel_pos,
        )
        for seq in sequences
        for mod in seq.modifications.get_mods(letter_mod_types)
        for match in seq.protein_matches[:None if all_matches else 1]
    )
