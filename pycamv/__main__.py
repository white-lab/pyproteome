
"""
Main module for running pycamv from the commandline.
"""

import argparse
import logging
import sys

from . import validate
from pyproteome import paths


def _parse_args(args):
    """
    Parses arguments.

    Parameters
    ----------
    args : list of str

    Returns
    -------
    argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser(
        prog="pycamv",
        description="Process searched mass spectrometry runs for CAMV.",
    )

    parser.add_argument(
        "basename",
        help="Base name of RAW and XML files.",
    )

    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Increase verbosity of output.",
    )

    parser.add_argument(
        "--base_dir",
        default=paths.BASE_DIR,
        help="Path to folder containing data files.",
    )

    parser.add_argument(
        "--scans",
        nargs="*", type=int,
        help="Individual scans to select for validation.",
    )

    return parser.parse_args(args)


def main(args):
    """
    Run spectra validation from commandline arguments.

    Parameters
    ----------
    args : list of str
    """
    args = _parse_args(args)

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)

    if args.base_dir:
        paths.set_base_dir(args.base_dir)

    validate.validate_spectra(args.basename, scan_list=args.scans)


if __name__ == "__main__":
    main(sys.argv[1:])
