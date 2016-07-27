"""
This module provides functionality for comparing spectra with a list of
fragments.
"""

from . import masses


CID_TOL = 1000
HCD_TOL = 10

COLLISION_TOLS = {
    "CID": CID_TOL,
    "HCD": HCD_TOL,
}


class PeptideHit:
    """
    Attributes
    ----------
    mz : float
    intensity : float
    name : str or None
    score : int or None
    predicted_mz : float or None
    match_list : dict of str, tuple of (float, float)
    num_losses : int
    """
    def __init__(
        self, mz, intensity,
        name=None, score=None, predicted_mz=None, match_list=None,
    ):
        self.mz = mz
        self.intensity = intensity
        self.name = name
        self.score = score
        self.predicted_mz = predicted_mz
        self.match_list = match_list
        self.num_losses = name.count("-") if name else 0

    @property
    def frac_error(self):
        if self.predicted_mz is None:
            return 0

        return abs(self.predicted_mz - self.mz) / self.predicted_mz


def _calculate_ion_score(ion_name):
    score = 0

    if ion_name.startswith("MH"):
        score += 12
    elif ion_name.startswith("b_") or ion_name.startswith("y_"):
        score += 10

    score -= ion_name.count("-")

    return score


def compare_spectra(
    spectra, frag_ions, charge, c13_num,
    tol=None,
):
    """
    Parameters
    ----------
    spectra : :class:`pymzml.spec.Spectrum<spec.Spectrum>`
    frag_ions : dict of str, float
    charge : int
    c13_num : int
    tol : float, optional

    Returns
    -------
    list of :class:`PeptideHit<pycamv.compare.PeptideHit>`
        Peak assignments for each peak in spectra.
    """
    if tol is None:
        tol = CID_TOL

    max_y = max(spectra.i)

    # C13 isotope calculations
    # XXX: Support C13 isotopes
    delta_c13 = masses.exact_mass({"C": [-1, 1]})
    delta_iso = [delta_c13/i for i in range(1, charge + 1)]

    # Iterate over each peak, looking for fragments to assign to it
    out = []

    # Reprofiled Peaks? Centroided Peaks?
    for mz, intensity in spectra.centroidedPeaks:
        # intensity, mz = np.array(intensity)
        peak_candidates = {
            ion_name: (ion_mz, abs(ion_mz - mz))
            for ion_name, ion_mz in frag_ions.items()
            if abs(ion_mz - mz) / ion_mz < 1.5 * tol / 1e6
        }

        if not peak_candidates:
            # XXX: AA Identification?
            out.append(PeptideHit(mz, intensity))
            continue

        # Calculate ion scores based on their names
        scores = {
            ion_name: _calculate_ion_score(ion_name)
            for ion_name in peak_candidates
        }

        # XXX: 0.5? 1
        closest_ion = min(peak_candidates, key=lambda x: peak_candidates[x][1])
        scores[closest_ion] += 1

        best_ion = max(scores, key=lambda x: scores[x])

        retain = (
            intensity > 0.1 * max_y or
            best_ion.startswith("MH") or
            best_ion.startswith("b_") or
            best_ion.startswith("y_") or
            (
                best_ion.count("-") == 0 and
                intensity > 0.025 * max_y
            )
        )

        if not retain:
            best_ion = None

        out.append(
            PeptideHit(
                mz,
                intensity,
                name=best_ion,
                score=scores.get(best_ion, None),
                predicted_mz=peak_candidates.get(best_ion, None),
                match_list=peak_candidates,
            )
        )

    # Free up memory so we don't keep duplicate info about a spectra's peaks
    del spectra._peaks
    del spectra._mz
    del spectra._i
    del spectra._centroidedPeaks
    del spectra._reprofiledPeaks

    return out
