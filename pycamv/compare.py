"""
This module provides functionality for comparing spectra with a list of
fragments.
"""

from . import fragments


CID_TOL = 1000
HCD_TOL = 10


class PeptideHit:
    def __init__(self, name, score, mz, predicted_mz):
        self.name = name
        self.score = score
        self.mz = mz
        self.predicted_mz = predicted_mz
        self.num_losses = name.count("-")


def compare_spectra(
    spectra, frag_ions, other_ions, charge, c13_num,
    tol=None,
):
    """
    Parameters
    ----------
    spectra : :class:`pymzml.spec.Spectrum`
    frag_ions : dict of str, float
    other_ions : dict of str, float
    charge : int
    c13_num : int
    tol : float, optional
    """
    if tol is None:
        tol = CID_TOL

    max_y = ...
    delta_iso = ...

    # XXX: AA Identification?
    # XXX: Support C13 isotopes
    pass
