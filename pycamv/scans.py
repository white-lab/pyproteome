"""
This module provides functionality for interacting with mass spec scan data.
"""

import logging
import re
import tempfile

from . import proteowizard


LOGGER = logging.getLogger("pycamv.scans")


class ScanQuery:
    """
    Attributes
    ----------
    scan : int
    isolation_mz : float or None
    window_offset : tuple of (int, int) or None
    precursor_scan : int or None
    collision_type : str or None
    """
    def __init__(
        self, scan,
        isolation_mz=None, window_offset=None, precursor_scan=None,
        collision_type=None,
    ):
        self.scan = scan
        self.precursor_scan = precursor_scan
        self.window_offset = window_offset
        self.isolation_mz = isolation_mz
        self.collision_type = collision_type


def _scanquery_from_spectrum(spectrum):
    """
    Parameters
    ----------
    spectrum : :class:`pymzml.scan.Spectrum<scan.Spectrum>`

    Returns
    -------
    :class:`ScanQuery<pycamv.validate.ScanQuery>`
    """
    prefix = {"mzml": "http://psi.hupo.org/ms/mzml"}

    scan = spectrum["id"]
    isolation_mz = spectrum["isolation window target m/z"]
    window_offset = (
        spectrum["isolation window lower offset"],
        spectrum["isolation window upper offset"],
    )
    collision_type = re.search(
        ".*@(\w+)\d+",
        spectrum["filter string"]
    ).groups(1).upper()

    spectrum_ref = spectrum.xmlTreeIterFree.find(
        "mzml:precursorList/mzml:precursor", prefix,
    ).get("spectrumRef")
    precursor_scan = re.search("scan=(\d+)", spectrum_ref).group(1)

    return ScanQuery(
        scan,
        precursor_scan=precursor_scan,
        window_offset=window_offset,
        isolation_mz=isolation_mz,
        collision_type=collision_type,
    )


def _c13_num(pep_query, scan_query):
    """
    Counts the number of C13 atoms in a query, based on the mass-error between
    the expected and isolated m/z values.

    Parameters
    ----------
    pep_query : :class:`PeptideQuery<pycamv.mascot.PeptideQuery>`
    scan_query : :class:`ScanQuery<pycamv.validate.ScanQuery>`

    Returns
    -------
    int
    """
    return int(
        round(
            pep_query.pep_exp_z *
            abs(pep_query.pep_exp_mz - scan_query.isolation_mz)
        )
    )


def get_scan_data(basename, pep_queries, out_dir=None):
    """
    Gets MS^2 and MS data for all scans in queries.

    Parameters
    ----------
    basename : str
    pep_queries : list of :class:`PeptideQuery<pycamv.mascot.PeptideQuery>`
    out_dir : str, optional

    Returns
    -------
    scan_queries : list of :class:`ScanQuery<pycamv.validate.ScanQuery>`
    ms2_data : :class:`pymzml.scan.Spectrum<scan.Spectrum>`
    ms_data : :class:`pymzml.scan.Spectrum<scan.Spectrum>`
    """
    if out_dir is None:
        out_dir = tempfile.mkdtemp()

    # Collect MS^2 data
    LOGGER.info("Converting MS^2 data.")
    ms2_data = proteowizard.raw_to_mzml(
        basename, out_dir,
        scans=sorted(set(pep_query.scan for pep_query in pep_queries)),
    )

    # Build a list of scan queries, including data about each scan
    # scan_queries = [
    #     _scanquery_from_spectrum(spectrum)
    #     for spectrum in ms2_data
    #     if spectrum["ms level"] == 2
    # ]
    scan_queries = [
        _scanquery_from_spectrum(ms2_data[pep_query.scan])
        for pep_query in pep_queries
    ]

    # Collect MS^1 data
    LOGGER.info("Converting MS^1 data.")
    ms_data = proteowizard.raw_to_mzml(
        basename, out_dir,
        scans=sorted(set(i.precursor_scan for i in scan_queries)),
    )

    return scan_queries, ms2_data, ms_data
