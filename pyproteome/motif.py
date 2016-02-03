"""
This module provides functionality for finding motifs in sequences.

Functionality includes n-mer generation.
"""

from __future__ import division

# Built-ins
import logging
import os
import random

# Core data analysis libraries
import Bio.Alphabet.IUPAC
import Bio.Seq
import Bio.motifs
import pandas as pd
from scipy.stats import hypergeom


LOGGER = logging.getLogger("pyproteome.sequence")


def _generate_n_mers(
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


class Motif:
    """
    Contains a motif that may match to one or more protein sequences.

    Attributes
    ----------
    motif : str

    Examples
    --------
    >>> import pyproteome
    >>> motif = pyproteome.Motif("O..x.-+")
    >>> "IEFyFER" in motif
    True
    >>> "IEFyFED" in motif
    False
    >>> "FFFFFFR" in motif
    False
    """
    def __init__(self, motif):
        self.motif = motif

    def generate_children(self):
        char_mapping = {
            "x": "sty",
            "O": "MILV",
            "-": "DE",
            "+": "RK",
            ".": "ACDEFGHIKLMNPQRSTVWYystO-+x",
        }
        for index, x in enumerate(self.motif):
            for char in char_mapping.get(x, ""):
                yield Motif(
                    self.motif[:index] + char + self.motif[index + 1:],
                )

    def __repr__(self):
        return "<pyproteome.Motif: {}>".format(self.motif)

    def __str__(self):
        return self.motif

    def __hash__(self):
        return hash(self.motif)

    def __contains__(self, other):
        if isinstance(other, Motif):
            other = other.motif

        return self._match(other)

    def __eq__(self, other):
        if isinstance(other, Motif):
            return self.motif == other.motif

        return self._match(other)

    matchers = {
        "x": lambda x: x in "sty",
        "O": lambda x: x in "MILV",
        "-": lambda x: x in "DE",
        "+": lambda x: x in "RK",
        ".": lambda x: True,
    }

    def _match(self, other):
        if not isinstance(other, str):
            raise TypeError

        if len(self.motif) != len(other):
            raise ValueError("Incompatible sequence lengths.")

        return all(
            Motif.matchers.get(motif, lambda x: x == motif)(char)
            for char, motif in zip(other, self.motif)
        )


def motif_enrichment(
    foreground, background,
    sig_cutoff=0.01, min_fore_hits=3, max_depth=None,
    start_letter="x", letter_mod_types=None,
    n_fpr_subsets=1000,
):
    """
    Calculate motifs significantly enriched in a list of sequences.

    Parameters
    ----------
    foreground : list of pyproteome.Sequence
    background : list of pyproteome.Sequence
    sig_cutoff : float, optional
    min_motifs : int, optional
    max_depth : int, optional
    start_letter : str, optional
    letter_mod_types : list of tuple of str, str

    Returns
    -------
    pandas.DataFrame

    Notes
    -----
    .. [1] Joughin, Brian a et al. "An Integrated Comparative Phosphoproteomic
           and Bioinformatic Approach Reveals a Novel Class of MPM-2 Motifs
           Upregulated in EGFRvIII-Expressing Glioblastoma Cells." Molecular
           bioSystems 5.1 (2009): 59-67.
    """
    if letter_mod_types is None:
        letter_mod_types = [(None, "Phospho")]

    foreground = list(
        _generate_n_mers(
            foreground,
            letter_mod_types=letter_mod_types,
        )
    )
    background = list(
        _generate_n_mers(
            background,
            letter_mod_types=letter_mod_types,
        )
    )
    assert len(foreground) > 0
    n = len(foreground[0])
    assert n % 2 == 1

    if max_depth is None:
        max_depth = n - 1

    LOGGER.info(
        "Starting analysis, n={}, N={}".format(
            len(foreground), len(background),
        )
    )
    LOGGER.info(
        "Motif sig cutoff: {}, Minimum hits: {}, Motif depth: {}".format(
            sig_cutoff, min_fore_hits, max_depth,
        )
    )

    def _motif_sig(fore_hits, back_hits):
        big_n = len(background)
        little_n = len(foreground)
        m = back_hits
        k = fore_hits
        return hypergeom.cdf(k, big_n, m, little_n)

    def _motif_stats(motif):
        fore_hits = sum(i in motif for i in foreground)
        back_hits = sum(i in motif for i in background)

        # Re-calculate the P enrichment score associated with this motif
        sig_p = _motif_sig(fore_hits, back_hits)

        # Calculate the p-value by comparing sig_p to the motif enrichment
        # score from random subsets of the background.
        p_value = sum(
            _motif_sig(
                sum(
                    i in motif
                    for i in random.sample(background, len(foreground))
                ),
                back_hits,
            ) <= sig_p
            for _ in range(n_fpr_subsets)
        ) / n_fpr_subsets

        return (motif, fore_hits, back_hits, sig_p, p_value)

    visited = set()
    motif_hits = set()

    # Set the starting motif and begin adding modifications to it.
    start = Motif(
        "." * (n // 2) + start_letter + "." * (n // 2)
    )

    LOGGER.info("Starting Motif: {}".format(start))

    to_process = [
        (motif, 1)
        for motif in start.generate_children()
    ]

    # First pass: Find all motifs with p < sig_cutoff
    for motif, depth in to_process:
        # Mark motif as visited
        visited.add(motif)

        # Calculate the number of foreground hits
        fore_hits = sum(i in motif for i in foreground)

        if fore_hits < min_fore_hits:
            continue

        # Calculate the number of background hits
        back_hits = sum(i in motif for i in background)

        # Calculate the motif signifiance
        sig_p = _motif_sig(fore_hits, back_hits)

        if sig_p >= sig_cutoff:
            continue

        # If we get this far, we have a hit!
        motif_hits.add(motif)

        # Now try adding to the list of modifications on the motif to find
        # something more specific.
        if depth >= max_depth:
            continue

        for new_motif in motif.generate_children():
            if new_motif in visited:
                continue

            to_process.append((new_motif, depth + 1))

    LOGGER.info("First pass: {} motifs".format(len(motif_hits)))

    # Second pass: Find all motifs that do not include a more specific subset
    # with the same number of foreground hits.
    motif_hits = set(
        motif
        for motif in motif_hits
        if not any(
            # other in motif: other is +more+ specific than motif
            other in motif and
            (
                sum(i in other for i in foreground) >=
                sum(i in motif for i in foreground)
            )
            for other in motif_hits
            if other != motif
        )
    )

    LOGGER.info(
        "Second pass (less specific filtered out): {} motifs"
        .format(len(motif_hits))
    )

    # Finally prepare the output as a sorted list with the motifs and their
    # associated statistics.
    LOGGER.info(
        "Calculating final p-values associated with motifs (N-FPR-subsets={})."
        .format(n_fpr_subsets)
    )
    motif_hits.add(start)

    df = pd.DataFrame(
        data=[_motif_stats(motif) for motif in motif_hits],
        columns=[
            "Motif",
            "Foreground Hits",
            "Background Hits",
            "p-score",
            "fpr-score",
        ],
    )
    df.sort(columns=["p-score"], inplace=True)
    df.reset_index(drop=True)

    return df


def motif_logo(
    data, letter_mod_types=None,
    folder_name=None, filename="Motif.svg",
):
    """
    Create a sequence logo figure.

    Logos are created based on the frequencies of peptides in a data set.

    Parameters
    ----------
    data : pyproteome.DataSet
    letter_mod_types : list of tuple of str, str
    folder_name : str, optional
    filename : str, optional
    """
    if folder_name and filename:
        filename = os.path.join(folder_name, filename)

    alpha = Bio.Alphabet.IUPAC.protein
    m = Bio.motifs.create(
        [
            Bio.Seq.Seq(seq.upper(), alphabet=alpha)
            for seq in _generate_n_mers(
                data.psms["Sequence"],
                letter_mod_types=letter_mod_types,
            )
        ],
        alphabet=alpha,
    )
    m.weblogo(filename, format="svg", logo_title="Testing")
