import os


def makedirs(folder_name=None):
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
