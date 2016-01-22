"""Utility functions used in other modules."""
# Built-ins
import difflib
from collections import OrderedDict
import os


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
        return [i + "_norm" for i in channels]

    if isinstance(channels, OrderedDict):
        return OrderedDict(
            (key + "_norm", val)
            for key, val in channels.items()
        )
