# -*- coding: utf-8 -*-

import argparse
import logging
import sys

from pycamv import validate, export


LOGGER = logging.getLogger("pycamv.main")

def _parse_args(args):
    """
    Parses arguments from a argv format.
    Parameters
    ----------
    args : list of str
    Returns
    -------
    argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser(
        prog="PyCAMV Converter",
        description="Aww yeah, mass specs!",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Increase verbosity of output.",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Turn on debugging mode.",
    )
    parser.add_argument(
        "--show_gui",
        action="store_true",
        help="Show GUI for converting files.",
    )
    parser.add_argument(
        "--basename",
    )
    parser.add_argument(
        "--raw_path",
    )
    parser.add_argument(
        "--xml_path",
    )
    parser.add_argument(
        "--scans_path",
    )
    parser.add_argument(
        "--out_path",
    )
    return parser.parse_args(args)


def run_gui():
    pass


def main(args):
    args = _parse_args(args)

    if args.verbose or args.debug:
        level = logging.DEBUG
    else:
        level = logging.INFO

    logging.basicConfig(
        filename="pycamv.log",
        level=level,
    )

    if args.show_gui:
        run_gui()
    else:
        if (
            args.basename is None or (
                args.xml_path is None and
                args.raw_path is None and
                args.out_dir is None
            )
        ):
            raise Exception(
                "Missing either basename or input xml / raw / scans / out path"
            )

        options, peak_hits, scan_mapping, precursor_windows, label_windows = (
            validate.validate_spectra(
                basename=args.basename,
                xml_path=args.xml_path,
                raw_path=args.raw_path,
                scans_path=args.scans_path,
            )
        )

        export.export_to_camv(
            args.out_path,
            peak_hits, scan_mapping, precursor_windows, label_windows,
        )

if __name__ == "__main__":
    try:
        main(sys.argv[1:])
    except Exception as e:
        LOGGER.error("PyCAMV Converter has crashed!", exc_info=True)
