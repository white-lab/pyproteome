'''Utility functions used in other modules.'''

# Built-ins
from collections import OrderedDict, Callable
import copy
import difflib
import functools
import os
import pickle
import types

import numpy as np
import pandas as pd

from . import paths


DEFAULT_DPI = 300
'''
The DPI to use when generating all image figures.
'''


def fuzzy_find(needle, haystack):
    '''
    Find the longest matching subsequence of needle within haystack.

    Returns the corresponding index from the beginning of needle.

    Parameters
    ----------
    needle : str
    haystack : str

    Returns
    -------
    index : int
    '''
    s = difflib.SequenceMatcher(a=haystack, b=needle)
    best = s.find_longest_match(0, len(haystack), 0, len(needle))
    return best.a - len(needle) + best.size


def make_folder(data=None, folder_name=None, sub='Output'):
    if folder_name is None:
        folder_name = os.path.join(
            paths.FIGURES_DIR,
            data.name
            if data is not None else
            'All',
            sub,
        )

    return makedirs(folder_name)


def makedirs(folder_name=None):
    '''
    Creates a folder if it does not exist.

    Parameters
    ----------
    folder_name : str, optional

    Returns
    -------
    folder_name : str
    '''
    if folder_name:
        try:
            os.makedirs(folder_name)
        except OSError:
            pass

    return folder_name


def norm(channels):
    '''
    Converts a list of channels to their normalized names.

    Parameters
    ----------
    channels : list of str or dict of (str, str) or None

    Returns
    -------
    new_channels : list of str or dict of str, str
    '''
    if channels is None:
        return None

    if isinstance(channels, str):
        return channels + '_norm'

    if isinstance(channels, list):
        return [norm(i) for i in channels]

    if isinstance(channels, (dict, OrderedDict)):
        return OrderedDict(
            (key, norm(val))
            for key, val in channels.items()
        )


def which(program):
    '''
    Checks if a program exists in PATH's list of directories.

    Parameters
    ----------
    program : str

    Returns
    -------
    path : str or None
    '''
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, _ = os.path.split(program)

    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ['PATH'].split(os.pathsep):
            path = path.strip('\'')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def flatten_set(lst):
    '''
    Flattens an Iterable with arbitrary nesting into a single set.

    Parameters
    ----------
    lst : Iterable

    Returns
    -------
    flattened : set

    Examples
    --------
        >>> utils.flatten_set([0, [1, 2], [[3]], 'string'])
        set([0, 1, 2, 3, 'string'])
    '''
    if isinstance(
        lst,
        (list, tuple, set, types.GeneratorType, pd.Series, np.ndarray)
    ):
        ret = set()

        for element in lst:
            for new_element in flatten_set(element):
                ret.add(new_element)

        return ret

    return set([lst])


def flatten_list(lst):
    '''
    Flattens an Iterable with arbitrary nesting into a single list.

    Parameters
    ----------
    lst : Iterable

    Returns
    -------
    flattened : list

    Examples
    --------
        >>> utils.flatten_list([0, [1, 2], [[3]], 'string'])
        [0, 1, 2, 3, 'string']
    '''
    if isinstance(
        lst,
        (list, tuple, set, types.GeneratorType, pd.Series, np.ndarray)
    ):
        ret = []

        for element in lst:
            for new_element in flatten_list(element):
                ret.append(new_element)

        return ret

    return [lst]


class DefaultOrderedDict(OrderedDict):
    # Source: http://stackoverflow.com/a/6190500/562769
    def __init__(self, default_factory=None, *a, **kw):
        if (default_factory is not None and
           not isinstance(default_factory, Callable)):
            raise TypeError('first argument must be callable')
        OrderedDict.__init__(self, *a, **kw)
        self.default_factory = default_factory

    def __getitem__(self, key):
        try:
            return OrderedDict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):
        if self.default_factory is None:
            args = tuple()
        else:
            args = self.default_factory,
        return type(self), args, None, None, self.items()

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        return type(self)(self.default_factory, self)

    def __deepcopy__(self, memo):
        return type(self)(self.default_factory,
                          copy.deepcopy(self.items()))

    def __repr__(self):
        return 'OrderedDefaultDict(%s, %s)' % (
            self.default_factory,
            OrderedDict.__repr__(self),
        )


def memoize(func):
    '''
    Memoize a function, saving its returned value for a given set of parameters
    in an in-memory cache.

    Examples
    --------
    >>> from pyproteome import utils
    >>> @utils.memoize
    ... def download_data(species):
    ...    ...  # Fetch / calculate the return value once


    Parameters
    ----------
    func : func

    Returns
    -------
    memorized : func
    '''
    cache = func.cache = {}

    @functools.wraps(func)
    def memoized_func(*args, **kwargs):
        key = str(args) + str(kwargs)

        if key not in cache:
            cache[key] = func(*args, **kwargs)

        return cache[key]

    return memoized_func


PICKLE_DIR = '.pyproteome'
'''
Default directory to use for saving / loading pickle files.
'''


def save(name, val=None):
    '''
    Save a variable using the pickle module.

    Parameters
    ----------
    name : str
        The name to use for data storage.
    val : object, optional

    Returns
    -------
    val : object
    '''
    filename = os.path.join(PICKLE_DIR, '{}.pkl'.format(name))

    makedirs(PICKLE_DIR)

    with open(filename, 'wb') as f:
        pickle.dump(val, f)

    return val


def load(name, default=None):
    '''
    Load a variable using the pickle module.

    Parameters
    ----------
    name : str
        The name to use for data storage.
    default : object, optional

    Returns
    -------
    val : object
    '''
    filename = os.path.join(PICKLE_DIR, '{}.pkl'.format(name))

    try:
        with open(filename, 'rb') as f:
            val = pickle.load(f)
    except (
        OSError, pickle.UnpicklingError, IOError,
        AttributeError, EOFError, ImportError, IndexError,
    ):
        val = default

    return val


def adjust_text(*args, **kwargs):
    '''
    Wraps importing and calling :func:`adjustText.adjust_text`.
    '''
    from adjustText import adjust_text as at
    return at(*args, **kwargs)


def get_name(proteins):
    '''
    Generates a shortened version of a protein name. For peptides
    that map to multiple proteins, this function finds the longest
    common prefix (excluding digits) that matches all proteins.

    Parameters
    ----------
    proteins : :class:`.data_sets.protein.Proteins`

    Returns
    -------
    str

    Examples
    --------
    >>> pyp.utils.get_name(
    ...     protein.Proteins([
    ...         protein.Protein(gene='Dpysl2'),
    ...         protein.Protein(gene='Dpysl3'),
    ...     ])
    ... )
    'Dpysl2/3'
    >>> pyp.utils.get_name(
    ...     protein.Proteins([
    ...         protein.Protein(gene='Src'),
    ...         protein.Protein(gene='Fgr'),
    ...         protein.Protein(gene='Fyn'),
    ...     ])
    ... )
    'Src / Fgr / Fyn'
    >>> pyp.utils.get_name(
    ...     protein.Proteins([
    ...         protein.Protein(gene='Tuba1a'),
    ...         protein.Protein(gene='Tuba1b'),
    ...         protein.Protein(gene='Tuba1c'),
    ...         protein.Protein(gene='Tuba4a'),
    ...         protein.Protein(gene='Tuba8'),
    ...     ])
    ... )
    'Tuba1a/1b/1c/3a/4a/8'
    '''
    genes = sorted(proteins.genes)
    common = ''
    sep = ' / '

    if len(genes) > 1:
        common = os.path.commonprefix(genes)
        last_digit = [
            ind
            for ind, i in list(enumerate(common))[::-1]
            if i.isdigit()
        ]
        
        if last_digit:
            common = common[:last_digit[0]]

        if common:
            sep = '/'

    return common + sep.join(i[len(common):] for i in genes)


def stars(p, ns='ns'):
    '''
    Calculate the stars to indicate significant changes.

    \\*\\*\\*\\* : p < 1e-4
    
    \\*\\*\\* : p < 1e-3
    
    \\*\\* : p < 1e-2
    
    \\* : p < 5e-2
    
    ns : not significant

    Parameters
    ----------
    p : float
    ns : str, optional

    Returns
    -------
    str
    '''
    if p < 1e-4:
        return '****'
    elif (p < 1e-3):
        return '***'
    elif (p < 1e-2):
        return '**'
    elif (p < 0.05):
        return '*'
    else:
        return ns
