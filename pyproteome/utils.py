"""Utility functions used in other modules."""

# Built-ins
from collections import OrderedDict
import difflib
import functools
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
            (norm(key), val)
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
