'''
This module provides functionality for manipulating sequences.
'''

# Built-ins
import logging

from . import modification, protein

import pyproteome as pyp


LOGGER = logging.getLogger('pyproteome.sequence')


class ProteinMatch:
    '''
    Contains information about how a peptide sequence maps onto a protein.

    Attributes
    ----------
    protein : :class:`.protein.Protein`
        Protein object.
    rel_pos : int
        Relative position of the peptide start within the protein sequence.
    exact : bool
        Indicates whether a peptide sequence exact matches its protein sequence.
    '''

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
    '''
    Contains information about a sequence and which proteins it matches to.

    Attributes
    ----------
    pep_seq : str
        Peptide sequence, in 1-letter amino code.
    protein_matches : list of :class:`.ProteinMatch`
        Object mapping all proteins that a peptide sequence matches.
    modifications : :class:`.modification.Modifications`
        Object listing all post-translation modifications identified on a peptide.
    '''

    def __init__(
        self,
        pep_seq='',
        protein_matches=None,
        modifications=None,
    ):
        '''
        Parameters
        ----------
        pep_seq : str
        protein_matches : list of :class:`.ProteinMatch`
        modifications : :class:`.modification.Modifications`, optional
        '''
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
            return other in [
                self.__str__(),
                self.__str__(show_mods=True),
                self.__str__(skip_labels=False, skip_terminus=False),
                self.__str__(
                    skip_labels=False,
                    skip_terminus=False,
                    show_mods=True,
                ),
            ]

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
            return any([
                other in i
                for i in [
                    self.__str__(),
                    self.__str__(skip_labels=True, skip_terminus=False),
                ]
            ])

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

    def __str__(
        self, 
        skip_labels=True, 
        skip_terminus=True, 
        mods=None, 
        show_mods=False,
    ):
        '''
        Converts a peptide into a string.

        Parameters
        ----------
        skip_labels : bool, optional
            Don't include TMT/iTRAQ quantitification tags in string.
        skip_terminus : bool, optional
            Don't show N-/C-terminal modifications.
        mods : list of str, optional
            Only show this subset of modifications (i.e. ['Phospho', 'Oxidation']).
        show_mods : bool, optional
            If true, show modification identities (i.e. 'y(Phospho)').
            Otherwise residues with modifications are shown as lowercase.

        Returns
        -------
        str
        '''
        string = list('N-' + self.pep_seq.upper() + '-C')
        self_mods = self.modifications

        if skip_terminus:
            string = string[2:-2]

        if mods:
            self_mods = self_mods.get_mods(mods)
        if skip_labels:
            self_mods = self_mods.skip_labels()

        def _mods(index, letter):
            lst = [
                mod
                for mod in self_mods
                if (mod.rel_pos == index)
                or (mod.nterm and index < 0)
                or (mod.cterm and index > len(self.pep_seq))
            ]

            if not lst or letter == '-':
                return letter

            if not show_mods:
                return letter.lower()

            return (
                letter.lower() +
                '({})'.format(', '.join([mod.mod_type for mod in lst]))
            )

        string = [
            _mods(ind, letter)
            for ind, letter in enumerate(
                string,
                start=0 if skip_terminus else -2,
            )
        ]
        return ''.join(string)

    def __len__(self):
        return len(self.pep_seq)

    @property
    def is_labeled(self):
        '''
        Checks whether a sequence is modified on any residue with a
        quantification label.

        Returns
        -------
        is_labeled : bool
        '''
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
        '''
        Checks whether a sequence is modified with quantification labels on
        fewer than all expected residues.

        Returns
        -------
        is_underlabeled : bool
        '''

        if self._is_underlabeled is not None:
            return self._is_underlabeled

        underlabeled = False

        if self.is_labeled:
            # XXX: Hardcodes label modification locations, not extendable to
            # new quantification tags without changes to this function
            underlabeled = not any(
                j.mod_type in modification.LABEL_NAMES and j.nterm
                for j in self.modifications.mods
            ) or self.pep_seq.count('K') != sum(
                j.mod_type in modification.LABEL_NAMES
                for j in self.modifications.mods
                if j.letter == 'K' and not j.nterm
            )

        self._is_underlabeled = underlabeled

        return underlabeled


def extract_sequence(proteins, sequence_string):
    '''
    Extract a Sequence object from a list of proteins and sequence string.

    Does not set the Sequence.modifications attribute.

    Parameters
    ----------
    proteins : list of :class:`.protein.Protein`
    sequence_string : str

    Returns
    -------
    seqs : list of :class:`.Sequence`
    '''
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
