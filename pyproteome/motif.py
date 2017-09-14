"""
This module provides functionality for finding motifs in sequences.

Functionality includes n-mer generation.
"""

from __future__ import division

# Built-ins
from collections import defaultdict
from functools import partial
import logging
import multiprocessing
import random
import re

# Core data analysis libraries
import pandas as pd
from scipy.stats import hypergeom

from . import sequence


LOGGER = logging.getLogger("pyproteome.motif")


class Motif:
    """
    Contains a motif that may match to one or more protein sequences.

    Matches include the regular single-letter amino acid names as well as
    phosphosites for serine, threonine, and tyrosine, non-polar amino acids,
    and positively and negatively charged amino acids.

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
        # Set _motif and _re explicitly, py2.7 does not call motif's setter
        # if we just set self.motif here...
        self._motif = motif
        self._re = self._compile_re(motif)

    char_mapping = {
        "x": "st",
        "O": "MILV",
        "-": "DE",
        "+": "RK",
        ".": "ACDEFGHIKLMNPQRSTVWYystO-+x",
    }

    @property
    def motif(self):
        return self._motif

    @motif.setter
    def motif(self, value):
        self._motif = value
        self._re = self._compile_re(value)

    def _compile_re(self, motif):
        ret = ""

        for char in motif:
            if char == "." or char not in self.char_mapping:
                ret += char
            else:
                ret += "[{}{}]".format(self.char_mapping[char], char)

        return re.compile(ret)

    def generate_children(self):
        for index, x in enumerate(self.motif):
            if x != ".":
                continue

            for char in self.char_mapping.get(x, ""):
                yield Motif(
                    self.motif[:index] + char + self.motif[index + 1:],
                )

    def generate_pairwise_children(self, hit_list):
        # First build a list of all of the AAs in the hits for this motif. We
        # don't want to generate child motifs that do not at least match any
        # peptides in that list.
        aa_is = dict(
            (
                (i, j),
                set((hit[i], hit[j]) for hit in hit_list[self]),
            )
            for i in range(len(self.motif))
            for j in range(i + 1, len(self.motif))
        )

        for key, val in aa_is.items():
            to_add = set()
            for special_char in "O-+x":
                for aa_i, aa_j in val:
                    if aa_i in self.char_mapping[special_char]:
                        to_add.add((special_char, aa_j))
                    if aa_j in self.char_mapping[special_char]:
                        to_add.add((aa_i, special_char))
                    if aa_i in self.char_mapping[special_char] and \
                       aa_j in self.char_mapping[special_char]:
                        to_add.add((special_char, special_char))
            val.update(to_add)

        # Then generate motifs with AAs that match the hit list
        for i, char_i in enumerate(self.motif):
            if char_i != ".":
                continue

            for j, char_j in enumerate(self.motif[i + 1:], start=i + 1):
                if char_j != ".":
                    continue

                for aa_i, aa_j in aa_is[(i, j)]:
                    yield Motif(
                        self.motif[:i] + aa_i +
                        self.motif[i + 1:j] + aa_j +
                        self.motif[j + 1:]
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

        if not isinstance(other, str):
            raise TypeError

        if len(self.motif) != len(other):
            raise ValueError("Incompatible sequence lengths.")

        return self.match(other)

    def __eq__(self, other):
        if isinstance(other, Motif):
            return self.motif == other.motif

        if not isinstance(other, str):
            raise TypeError

        if len(self.motif) != len(other):
            raise ValueError("Incompatible sequence lengths.")

        return self.match(other)

    def __lt__(self, other):
        if not isinstance(other, Motif):
            raise TypeError

        return self.motif < other.motif

    def match(self, other):
        """
        Match a given sequence to this motif.

        Parameters
        ----------
        other : str

        Returns
        -------
        bool
        """
        return bool(self._re.match(other))


def get_nmer_args(kwargs):
    nmer_args = {}

    nmer_args["letter_mod_types"] = kwargs.pop(
        "letter_mod_types", [(None, "Phospho")],
    )
    nmer_args["use_ptms"] = kwargs.pop("use_ptms", True)
    nmer_args["use_nterms"] = kwargs.pop("use_nterms", False)
    nmer_args["use_cterms"] = kwargs.pop("use_cterms", False)

    return nmer_args


def generate_n_mers(
    sequences, n=15,
    all_matches=True,
    fill_left="A",
    fill_right="A",
    letter_mod_types=None,
    use_ptms=True,
    use_nterms=False,
    use_cterms=False,
):
    """
    Generate n-mers around all sites of modification in sequences.

    Parameters
    ----------
    sequences : list of :class:`Sequence<pyproteome.sequence.Sequence>`
    n : int, optional
    all_matches : bool, optional
        Generate n-mers for all protein matches else just the first match.
    fill_left : str, optional
    fill_right : str, optional
    letter_mod_types : list of tuple of str, str, optional
    use_ptms : bool, optional
    use_nterms : bool, optional
    use_cterms : bool, optional
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

    def _get_seqs(seq):
        for match in seq.protein_matches[:None if all_matches else 1]:
            if use_ptms:
                for mod in seq.modifications.get_mods(letter_mod_types):
                    yield match.rel_pos + mod.rel_pos, match

            if use_nterms:
                yield match.rel_pos, match

            if use_cterms:
                yield match.rel_pos + len(seq.pep_seq) - 1, match

    return set(
        _n_mer_from_sequence(
            match.protein.full_sequence,
            pos,
        )
        for seq in sequences
        for pos, match in _get_seqs(seq)
    )


def _motif_stats(motif, fore_size, back_size, done, pp_dist=None):
    fore_hits, back_hits = done[motif]

    # Re-calculate the P enrichment score associated with this motif
    p_value = _motif_sig(fore_hits, fore_size, back_hits, back_size)

    ppvalue = _pp_sig(p_value, pp_dist) if pp_dist else None

    return (
        motif,
        fore_hits, fore_size,
        back_hits, back_size,
        p_value, ppvalue,
    )


def _is_a_submotif_with_same_size(submotif, fore_hits, done):
    """
    Test whether the 'submotif' is an ancestor of any known motif, and
    matches the same number of foreground sequences.
    """
    # We want to know if this "submotif" contains any completed motifs
    return any(
        checked in submotif
        for checked, checked_hits in done.items()
        if fore_hits == checked_hits[0] and checked != submotif
    )


def _make_nmers(lst, motif_length, letter_mod_types):
    if isinstance(next(iter(lst)), sequence.Sequence):
        return list(
            generate_n_mers(
                lst,
                n=motif_length,
                letter_mod_types=letter_mod_types,
            )
        )

    assert all(len(i) == motif_length for i in lst)
    return lst


def _filter_less_specific(prev_pass, done):
    return set(
        motif
        for motif in prev_pass
        if not any(
            # other in motif: other is +more+ specific than motif
            other in motif and done[other][0] >= done[motif][0]
            for other in prev_pass
            if other != motif
        )
    )


def _check_anchestry(motif, parent, done):
    motif_length = len(motif.motif)

    for index, char in enumerate(motif.motif):
        if index == motif_length // 2 or char != ".":
            continue

        new_motif = Motif(
            motif.motif[:index] + "." + motif.motif[index + 1:]
        )

        if new_motif != parent and new_motif in done:
            return True

    return False


def _count_occurences(motif, parent, hit_list):
    hits = [
        i
        for i in hit_list[parent]
        if motif.match(i)
    ]
    hit_list[motif] = hits

    return len(hits)


def _motif_sig(fore_hits, fore_size, back_hits, back_size):
    return (
        1 -
        hypergeom.cdf(fore_hits, back_size, back_hits, fore_size) +
        hypergeom.pmf(fore_hits, back_size, back_hits, fore_size)
    )


def _pp_sig(p_value, pp_dist):
    return len([i for i in pp_dist if i < p_value]) / len(pp_dist)


# Taken from https://stackoverflow.com/questions/1023038/
def lowpriority():
    """ Set the priority of the process to below-normal."""

    import sys
    try:
        sys.getwindowsversion()
    except AttributeError:
        isWindows = False
    else:
        isWindows = True

    if isWindows:
        # Based on:
        #   "Recipe 496767: Set Process Priority In Windows" on ActiveState
        #   http://code.activestate.com/recipes/496767/
        import win32api
        import win32process
        import win32con

        pid = win32api.GetCurrentProcessId()
        handle = win32api.OpenProcess(win32con.PROCESS_ALL_ACCESS, True, pid)
        win32process.SetPriorityClass(
            handle, win32process.BELOW_NORMAL_PRIORITY_CLASS,
        )
    else:
        import os

        os.nice(1)


def _random_pdist(x, background, fore_size, kwargs):
    lowpriority()

    return motif_enrichment(
        random.sample(background, fore_size),
        background,
        **kwargs
    )[1]


def _generate_ppdist(
    background, fore_size, pp_iterations,
    cpu_count=None,
    **kwargs
):
    if cpu_count is None:
        try:
            cpu_count = multiprocessing.cpu_count()
        except NotImplementedError:
            cpu_count = 2

    kwargs = kwargs.copy()
    kwargs["pp_value"] = False

    pp_dist = []

    pool = multiprocessing.Pool(
        processes=cpu_count,
    )

    for p_dist in pool.imap_unordered(
        partial(
            _random_pdist,
            background=background,
            fore_size=fore_size,
            kwargs=kwargs,
        ),
        range(pp_iterations),
    ):
        pp_dist.append(
            min(p_dist + [kwargs.get("sig_cutoff")])
        )

    return pp_dist


def motif_enrichment(
    foreground, background,
    sig_cutoff=0.01, min_fore_hits=0,
    start_letters=None, letter_mod_types=None,
    motif_length=15,
    pp_value=False, pp_iterations=100,
    cpu_count=None,
):
    """
    Calculate motifs significantly enriched in a list of sequences.

    Parameters
    ----------
    foreground : list of :class:`Sequence<pyproteome.sequence.Sequence>`
    background : list of :class:`Sequence<pyproteome.sequence.Sequence>`
    sig_cutoff : float, optional
    min_fore_hits : int, optional
    start_letters : list of str, optional
    letter_mod_types : list of tuple of str, str
    motif_length : int, optional
    pp_value : bool, optional
    pp_iterations : int, optional
    cpu_count : int, optional
        Number of CPUs to use when calculating pp-values, does not apply to
        a single motif-enrichment process.

    Returns
    -------
    df : :class:`pandas.DataFrame`
    p_dist : list of float
    pp_dist : list of float

    Notes
    -----
    .. [1] Joughin, Brian a et al. "An Integrated Comparative Phosphoproteomic
           and Bioinformatic Approach Reveals a Novel Class of MPM-2 Motifs
           Upregulated in EGFRvIII-Expressing Glioblastoma Cells." Molecular
           bioSystems 5.1 (2009): 59-67.
    """
    p_dist = []

    def _search_children(children, parent=None):
        """
        Run a depth-first search over a given motif.
        """
        motif_hits = set()

        for motif in children:
            if motif in visited:
                failed["checkpat"] += 1
                continue

            # Mark motif as visited
            visited.add(motif)

            if motif in done:
                failed["done"] += 1
                continue

            # Calculate the number of foreground hits
            fore_hits = _count_occurences(motif, parent, fg_hit_list)

            if fore_hits < min_fore_hits:
                failed["count"] = 1
                del fg_hit_list[motif]
                continue

            # Check if we've already completed the search for a parent of this
            # motif
            if parent and _check_anchestry(motif, parent, done):
                failed["anchestry"] += 1
                del fg_hit_list[motif]
                continue

            # Shortcut calculating back-hits if we can help it
            best_p = _motif_sig(fore_hits, fore_size, fore_hits, back_size)
            if best_p > sig_cutoff:
                failed["bestsig"] += 1
                del fg_hit_list[motif]
                continue

            if _is_a_submotif_with_same_size(motif, fore_hits, done):
                failed["sub"] += 1
                del fg_hit_list[motif]
                continue

            # Calculate the number of background hits
            back_hits = _count_occurences(motif, parent, bg_hit_list)

            # Check the signifiance of the motif
            p_value = _motif_sig(fore_hits, fore_size, back_hits, back_size)
            p_dist.append(p_value)

            if p_value > sig_cutoff:
                failed["sig"] += 1
                del fg_hit_list[motif]
                del bg_hit_list[motif]
                continue

            # If we get this far, we have a hit!
            motif_hits.add(motif)
            failed["succeed"] += 1
            done[motif] = (fore_hits, back_hits)

            # Now try adding to the list of modifications on the motif to find
            # something more specific.
            motif_hits.update(
                _search_children(
                    children=motif.generate_children(),
                    parent=motif,
                )
            )
            motif_hits.update(
                _search_children(
                    children=motif.generate_pairwise_children(fg_hit_list),
                    parent=motif,
                )
            )

            del fg_hit_list[motif]
            del bg_hit_list[motif]

        return motif_hits

    fore_size = len(foreground)
    back_size = len(background)

    assert fore_size > 0
    assert back_size > 0
    assert motif_length % 2 == 1

    if letter_mod_types is None:
        letter_mod_types = [(None, "Phospho")]

    foreground = _make_nmers(foreground, motif_length, letter_mod_types)
    background = _make_nmers(background, motif_length, letter_mod_types)
    fore_size = len(foreground)
    back_size = len(background)

    LOGGER.info(
        "Starting analysis, n={}, N={}".format(
            fore_size, back_size,
        )
    )
    LOGGER.info(
        "Motif sig cutoff: {}, Minimum hits: {}".format(
            sig_cutoff, min_fore_hits,
        )
    )

    visited, done = set(), {}
    fg_hit_list = defaultdict(list)
    bg_hit_list = defaultdict(list)
    failed = defaultdict(int)

    # Set the starting motif and begin adding modifications to it.
    # Note: Order matters, put less specific to the end of this list
    if start_letters is None:
        start_letters = ["k", "y", "s", "t", "x"]

    starts = [
        Motif(
            "." * (motif_length // 2) + letter + "." * (motif_length // 2)
        )
        for letter in start_letters
    ]

    LOGGER.info(
        "Starting Motifs: {}".format(", ".join(str(i) for i in starts))
    )

    for start in starts:
        fg_hit_list[start] = foreground
        bg_hit_list[start] = background

    first_pass = set()

    try:
        for start in starts:
            first_pass.update(
                _search_children(
                    children=start.generate_children(),
                    parent=start,
                )
            )

            first_pass.update(
                _search_children(
                    children=start.generate_pairwise_children(fg_hit_list),
                    parent=start,
                )
            )

            del fg_hit_list[start]
            del bg_hit_list[start]
    except KeyboardInterrupt:
        return
    finally:
        # Output some debugging information to compare with Brian's output
        LOGGER.info(
            "\n".join(
                "FAIL_{}: {}".format(key.upper(), failed[key])
                for key in [
                    "done", "sub", "count", "sig", "bestsig",
                    "anchestry", "checkpat",
                ]
            ) + "\nSUCCEED: {}".format(failed["succeed"])
        )

    LOGGER.info("First pass: {} motifs".format(len(first_pass)))

    second_pass = _filter_less_specific(first_pass, done)

    LOGGER.info(
        "Second pass (less specific motifs removed): {} motifs"
        .format(len(second_pass))
    )

    pp_dist = _generate_ppdist(
        background, fore_size, pp_iterations,
        sig_cutoff=sig_cutoff, min_fore_hits=min_fore_hits,
        start_letters=start_letters, letter_mod_types=letter_mod_types,
        motif_length=motif_length,
        cpu_count=cpu_count,
    ) if pp_value else None

    # Finally prepare the output as a sorted list with the motifs and their
    # associated statistics.
    df = pd.DataFrame(
        data=[
            _motif_stats(motif, fore_size, back_size, done, pp_dist=pp_dist)
            for motif in second_pass
        ],
        columns=[
            "Motif",
            "Foreground Hits",
            "Foreground Size",
            "Background Hits",
            "Background Size",
            "p-value",
            "pp-value",
        ],
    )
    df.sort_values(by=["pp-value", "p-value", "Motif"], inplace=True)
    df.reset_index(drop=True, inplace=True)

    return df, p_dist, pp_dist
