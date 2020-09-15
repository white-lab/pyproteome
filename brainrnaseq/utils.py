
import functools
import os


def makedirs(folder_name=None):
    '''
    Creates a folder if it does not exist.

    Parameters
    ----------
    folder_name: str or None
    '''
    if folder_name:
        try:
            os.makedirs(folder_name)
        except OSError:
            pass


def memoize(func):
    cache = func.cache = {}

    @functools.wraps(func)
    def memoized_func(*args, **kwargs):
        key = str(args) + str(kwargs)

        if key not in cache:
            cache[key] = func(*args, **kwargs)

        return cache[key]

    return memoized_func
