"""
This module provides functionality for interacting with CAMV.

Currently limited to importing and outputing scan lists.
"""

# Built-ins
from __future__ import division

import logging
import re
import tempfile
import shutil

# Core data analysis libraries

from . import mascot, ms_labels, proteowizard


LOGGER = logging.getLogger("pycamv.validate")


class SearchOptions:
    """
    Contains options used to search a data set in MASCOT.

    Attributes
    ----------
    label_type : tuple of (str, int)
    """
    def __init__(self, fixed_mods, var_mods):
        self.label_type = (None, 0)

        for mod in fixed_mods:
            mod = re.sub(r"([\w-]+) \([\w-]+\)", r"\1", mod)
            num = ms_labels.LABEL_NUMBERS.get(mod, 0)

            if num > 0:
                self.label_type = ("Fixed", num)

        for mod in var_mods:
            mod = re.sub(r"([\w-]+) \([\w-]+\)", r"\1", mod)
            num = ms_labels.LABEL_NUMBERS.get(mod, 0)

            if num > 0:
                self.label_type = ("Variable", num)

        # TODO: Parse out SILAC, C-mod, phospho, etc


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


def validate_spectra(basename, scan_list=None):
    """
    Generate CAMV web page for validating spectra.

    Parameters
    ----------
    basename : str
    scan_list : list of int, optional
    """
    # Read MASCOT xml file
    xml_name = "{}.xml".format(basename)

    fixed_mods, var_mods, pep_queries = mascot.read_mascot_xml(xml_name)

    # Optionally filter queries using a scan list
    if scan_list is not None:
        pep_queries = [
            query
            for query in pep_queries
            if query.scan in scan_list
        ]

    # Extract MASCOT search options
    options = SearchOptions(fixed_mods, var_mods)

    # Remove peptides with an excess of modification combinations

    # Get scan data from RAW file
    out_dir = tempfile.mkdtemp()

    scan_queries, ms2_data, ms_data = get_scan_data(
        basename, pep_queries, out_dir,
    )

    # Determine SILAC precursor masses

    # Get Precursor Scan information

    # Remove precursor contaminated scans from validation list

    # Check for Cysteine carbamidomethylation present in MASCOT search

    # Check each assignment to each scan

    # Output data

    del ms_data
    del ms2_data
    shutil.rmtree(out_dir)

    return options


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
    list of :class:`ScanQuery<pycamv.validate.ScanQuery>`
    """
    if out_dir is None:
        out_dir = tempfile.mkdtemp()

    # Collect MS^2 data
    LOGGER.info("Converting MS^2 data.")
    ms2_data = proteowizard.raw_to_mzml(
        basename, out_dir,
        scans=[pep_query.scan for pep_query in pep_queries],
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
