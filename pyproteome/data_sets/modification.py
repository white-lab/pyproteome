'''
This module provides functionality for post-translational modifications.

Wraps modifications in a structured class and allows filtering of
modifications by amino acid and modification type.
'''

from collections import defaultdict
import copy

LABEL_NAME_TARGETS = (
    'TMT', 'ITRAQ', 'plex',
)
'''
Substrings used to identify and import novel label names from .msf files.
'''

LABEL_NAMES = defaultdict(set)
LABEL_NAMES['TMT10'].add('K')
LABEL_NAMES['TMT10'].add('N-term')
LABEL_NAMES['TMT6'].add('K')
LABEL_NAMES['TMT6'].add('N-term')
'''
Names of modifications used for quantification of peptide abundances.
'''
MERGE_UNDERLABELED = True
'''
Merge peptides that have satured TMT labeling with peptides that are underlabeled.
'''


class Modifications:
    '''
    A list of modifications.

    Wraps the Modification objects and provides several utility functions.

    Attributes
    ----------
    mods : list of :class:`.Modification`
    '''

    def __init__(self, mods=None):
        '''
        Initialize from a list of modifications.

        Parameters
        ----------
        mods : list of :class:`.Modification`
        '''
        self.mods = mods or ()

    def __iter__(self):
        return iter(self.mods)

    def __len__(self):
        return len(self.mods)

    def copy(self):
        '''
        Creates a copy of a set of modifications. Does not copy the underlying
        sequence object.

        Returns
        -------
        mods : :class:`.Modifications`
        '''
        new = copy.copy(self)
        new.mods = tuple(i.copy() for i in new.mods)
        return new

    def skip_labels(self):
        '''
        Get modifications, skipping over any that are peptide labels.

        Returns
        -------
        mods : list of :class:`.Modification`
        '''
        return [
            mod
            for mod in self.mods
            if not any(label in mod.mod_type for label in LABEL_NAMES)
        ]

    def get_mods(self, letter_mod_types):
        '''
        Filter the list of modifications.

        Only keeps modifications with a given letter, mod_type, or both.

        Examples
        --------
        >>> from pyproteome.sequence import Sequence
        >>> from pyproteome.modification import Modification, Modifications
        >>> s = Sequence(pep_seq='SVYTEIK')
        >>> m = Modifications(
        ...     [
        ...         Modification(mod_type='TMT', nterm=True, sequence=s),
        ...         Modification(mod_type='Phospho', rel_pos=2, sequence=s),
        ...         Modification(mod_type='TMT', rel_pos=6, sequence=s),
        ...     ]
        ... )
        >>> m.get_mods('TMT')
        ['TMT A0', 'TMT K6']
        >>> m.get_mods('Phospho')
        ['pY2']
        >>> m.get_mods('Y')
        ['pY2']
        >>> m.get_mods('S')
        []
        >>> m.get_mods([('Y', 'Phospho')])
        ['pY2']
        >>> m.get_mods([('S', 'Phospho')])
        []

        Parameters
        ----------
        letter_mod_types : list of tuple of str, str

        Returns
        -------
        mods : :class:`.Modifications`
        '''
        any_letter, any_mod, letter_mod = \
            _extract_letter_mods(letter_mod_types)
        return Modifications(
            tuple(
                mod
                for mod in self.mods
                if allowed_mod_type(
                    mod,
                    any_letter=any_letter,
                    any_mod=any_mod,
                    letter_mod=letter_mod,
                )
            )
        )

    def __hash__(self):
        return hash(
            tuple(
                sorted(self.skip_labels(), key=lambda x: x.to_tuple())
                # sorted(self.mods, key=lambda x: x.to_tuple())
            ),
        )

    def __eq__(self, other):
        if not isinstance(other, Modifications):
            raise TypeError()

        if MERGE_UNDERLABELED:
            self_mods = sorted(self.skip_labels(), key=lambda x: x.to_tuple())
            o_mods = sorted(other.skip_labels(), key=lambda x: x.to_tuple())
        else:
            self_mods = sorted(self.mods, key=lambda x: x.to_tuple())
            o_mods = sorted(other.mods, key=lambda x: x.to_tuple())

        return tuple(self_mods) == tuple(o_mods)

    def __lt__(self, other):
        return self.skip_labels() < other.skip_labels()

    def __repr__(self, absolute=True, skip_labels=True):
        return self.__str__(absolute=absolute, skip_labels=skip_labels)

    def __str__(
        self, 
        absolute=True, 
        skip_labels=True, 
        prot_index=None, 
        show_mod_type=True,
    ):
        if len(self.mods) == 0:
            return ''

        if skip_labels:
            lst = list(self.skip_labels())
        else:
            lst = list(iter(self))

        if not lst:
            return ''

        def _mod_prot(i):
            return ', '.join(
                '{}{}{}{}'.format(
                    mod.display_mod_type() if show_mod_type else '',
                    mod.letter,
                    1 + (mod.abs_pos[i] if absolute else mod.rel_pos),
                    '' if mod.exact[i] else '*'
                )
                for mod in lst
            )

        if prot_index is None:
            return ' / '.join(
                _mod_prot(i)
                for i in range(len(lst[0].exact))
            )
        else:
            return _mod_prot(prot_index)


class Modification:
    '''
    Contains information for a single peptide modification.

    Attributes
    ----------
    rel_pos : int
        The relative position of a modification in a peptide sequence
        (0-indexed).
    mod_type : str
        A short name for this type of modification (i.e. 'Phospho',
        'Carbamidomethyl', 'Oxidation', 'TMT6', 'TMT10')
    nterm : bool
        Boolean indicator of whether this modification is applied to the
        peptide N-terminus.
    cterm : bool
        Boolean indicator of whether this modification is applied to the
        peptide C-terminus.
    '''

    def __init__(
        self,
        rel_pos=0,
        mod_type='',
        sequence=None,
        nterm=False,
        cterm=False,
    ):
        self.rel_pos = rel_pos
        self.mod_type = mod_type
        self.nterm = nterm
        self.cterm = cterm
        self.sequence = sequence

    def display_mod_type(self):
        '''
        Return the mod_type in an abbreviated form (i.e. 'p' for 'Phospho')

        Returns
        -------
        abbrev : str
        '''
        if self.mod_type in ['Phospho']:
            return 'p'
        if self.mod_type in ['Carbamidomethyl']:
            return 'cm'
        if self.mod_type in ['Oxidation']:
            return 'ox'

        return self.mod_type

    def to_tuple(self):
        return (
            self.rel_pos,
            self.mod_type,
            self.nterm,
            self.cterm,
            self.letter,
            self.abs_pos,
            self.exact,
        )

    def __hash__(self):
        return hash(self.to_tuple())

    def __lt__(self, other):
        return self.to_tuple() < other.to_tuple()

    def __eq__(self, other):
        if not isinstance(other, Modification):
            raise TypeError()

        return self.to_tuple() == other.to_tuple()

    def copy(self):
        '''
        Creates a copy of a modification. Does not copy the underlying sequence
        object.

        Returns
        -------
        mod : :class:`.Modification`
        '''
        new = copy.copy(self)
        return new

    @property
    def letter(self):
        '''
        This modification's one-letter amino acid code (i.e. 'Y'), or 'N-term'
        / 'C-term' for terminal modifications.

        Returns
        -------
        letter : str
        '''
        if self.sequence is None:
            return ''

        if self.nterm:
            return 'N-term'
        elif self.cterm:
            return 'C-term'

        return self.sequence.pep_seq[self.rel_pos].upper()

    @property
    def abs_pos(self):
        '''
        The absolute positions of this modification in the full sequence
        of each mapped protein (0-indexed).

        Returns
        -------
        tuple of int
        '''
        if self.sequence is None:
            return ()

        return tuple(
            self.rel_pos + match.rel_pos
            for match in self.sequence.protein_matches
        )

    @property
    def exact(self):
        '''
        Indicates whether each peptide-protein mapping for this modification is
        an exact or partial match.

        Returns
        -------
        exact : tuple of bool
        '''
        if self.sequence is None:
            return ()

        return tuple(
            match.exact
            for match in self.sequence.protein_matches
        )

    def __repr__(self):
        return (
            '<Modification {}{}({})>'
        ).format(
            self.letter,
            (self.rel_pos + 1) if not self.cterm and not self.nterm else '',
            self.mod_type,
        )


def allowed_mod_type(mod, any_letter=None, any_mod=None, letter_mod=None):
    '''
    Check if a modification is of a given type.

    Filters by letter, mod_type, or both.

    Parameters
    ----------
    mod : :class:`.Modification`
    any_letter : set of str
    any_mod : set of str
    letter_mod : set of tuple of str, str

    Returns
    -------
    is_type : bool
    '''
    return (
        (
            any_letter is None or
            mod.letter in any_letter
        ) or (
            any_mod is None or
            mod.mod_type in any_mod
        ) or (
            letter_mod is None or
            (mod.letter, mod.mod_type) in letter_mod
        )
    )


def _extract_letter_mods(letter_mod_types=None):
    if letter_mod_types is None:
        return None, None, None

    if isinstance(letter_mod_types, str):
        letter_mod_types = (letter_mod_types,)

    any_letter = set()
    any_mod = set()
    letter_mod = set()

    for elem in letter_mod_types:
        if not isinstance(elem, tuple):
            if len(elem) == 1:
                any_letter.add(elem.upper())
            else:
                any_mod.add(elem)

            continue

        letter, mod_type = elem

        if letter is None and mod_type is None:
            raise Exception('Need at least one letter or mod type not None')
        elif letter is None and mod_type is not None:
            any_mod.add(mod_type)
        elif letter is not None and mod_type is None:
            any_letter.add(letter.upper())
        else:
            letter_mod.add((letter.upper(), mod_type))

    return any_letter, any_mod, letter_mod
