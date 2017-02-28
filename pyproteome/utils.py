"""Utility functions used in other modules."""

# Built-ins
from collections import OrderedDict, Callable
import copy
import difflib
import os

import numpy as np
import pandas as pd


def fuzzy_find(needle, haystack):
    """
    Find the longest matching subsequence of needle within haystack.

    Returns the corresponding index from the beginning of needle.

    Parameters
    ----------
    needle : str
    haystack : str

    Returns
    -------
    int
    """
    s = difflib.SequenceMatcher(a=haystack, b=needle)
    best = s.find_longest_match(0, len(haystack), 0, len(needle))
    return best.a - len(needle) + best.size


def make_folder(folder_name=None):
    """
    Creates a folder if it does not exist.

    Parameters
    ----------
    folder_name: str or None
    """
    if folder_name:
        try:
            os.makedirs(folder_name)
        except OSError:
            pass


def norm(channels):
    """
    Converts a list of channels to their normalized names.

    Parameters
    ----------
    channels : list of str or OrderedDict of str, str or None

    Returns
    -------
    list of str or OrderedDict of str, str
    """
    if channels is None:
        return None

    if isinstance(channels, str):
        return channels + "_norm"

    if isinstance(channels, list):
        return [norm(i) for i in channels]

    if isinstance(channels, OrderedDict):
        return OrderedDict(
            (key, norm(val))
            for key, val in channels.items()
        )


def which(program):
    """
    Checks if a program exists in PATH's list of directories.

    Parameters
    ----------
    program : str

    Returns
    -------
    str or None
    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)

    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def flatten_set(lst):
    if isinstance(lst, (list, tuple, set, pd.Series, np.ndarray)):
        ret = set()

        for element in lst:
            for new_element in flatten_set(element):
                ret.add(new_element)

        return ret

    return set([lst])


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
        return 'OrderedDefaultDict(%s, %s)' % (self.default_factory,
                                               OrderedDict.__repr__(self))
