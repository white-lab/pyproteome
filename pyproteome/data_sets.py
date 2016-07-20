"""
This module provides functionality for manipulating proteomics data sets.

Functionality includes merging data sets and interfacing with attributes in a
structured format.
"""

# Built-ins
from collections import OrderedDict
import copy
import logging

# Core data analysis libraries
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind

from . import loading, utils


LOGGER = logging.getLogger("pyproteome.data_sets")


class DataSet:
    """
    Class that encompasses a proteomics data set.

    Includes peptide list, scan list, channels, groups, quantitative
    phenotypes, etc.

    Attributes
    ----------
    psms : :class:`pandas.DataFrame`
        Contains at least "Proteins", "Sequence", and "Modifications" columns.
    channels : dict of str, str
        Maps label channel to sample name.
    groups : dict of str, list of str
        Maps groups to list of sample names. The primary group is considered
        as the first in this sequence.
    phenotypes : dict of str, (dict of str, float)
        Primary key is the phenotype and its value is a dictionary mapping
        sample names to that phenotype's value.
    name : str
    enrichments : list of str
    tissues : list of str
    sets : int
        Number of sets merged into this data set.
    source : str
    scan_list : dict of str, list of int
    """
    def __init__(
        self, channels,
        psms=None,
        mascot_name=None, camv_name=None, msf=True,
        groups=None, phenotypes=None,
        name="", enrichments=None, tissues=None,
        dropna=True,
        camv_slices=None,
        merge_duplicates=True,
        merge_subsets=False,
        filter_bad=True,
    ):
        """
        Initializes a data set.

        Parameters
        ----------
        channels : dict of str, str
        psms : :class:`pandas.DataFrame`, optional
            Read psms directly from a DataFrame object.
        mascot_name : str, optional
            Read psms from MASCOT / Discoverer data files.
        camv_name : str, optional
            Read psms from CAMV data files.
        msf : bool, optional
            Read mascot data from .msf file instead of a tab-delimited file.
        groups : dict of str, list of str, optional
        phenotypes : dict of str, (dict of str, float), optional
        name : str, optional
        enrichments : list of str, optional
        tissues : list of str, optional
        dropna : bool, optional
        camv_slices : int, optional
        merge_duplicates : bool, optional
        merge_subsets : bool, optional
        filter_bad : bool, optional
        """
        assert (
            psms is not None or
            mascot_name is not None or
            camv_name is not None
        )

        self.source = "unknown"
        self.scan_lists = None
        self.validated = False

        if mascot_name:
            psms, self.scan_lists = loading.load_mascot_psms(
                mascot_name,
                camv_slices=camv_slices,
                msf=msf,
            )
            self.source = "MASCOT"
        elif camv_name:
            psms = loading.load_validated_psms(
                camv_name,
            )
            self.source = "CAMV"

        self.psms = psms
        self.channels = channels
        self.groups = groups
        self.phenotypes = phenotypes
        self.name = name
        self.enrichments = enrichments
        self.tissues = tissues

        self.normalized = False
        self.sets = 1

        if dropna:
            LOGGER.info("Dropping channels with NaN values.")
            self.dropna(inplace=True)

        if filter_bad:
            # self.psms = self.psms[
            #     self.psms["Confidence Level"].isin(["Medium", "High"])
            # ]
            self.psms = self.psms[
                self.psms["Confidence Level"].isin(["High"])
            ]

        if merge_duplicates:
            LOGGER.info("Merging duplicate peptide hits together.")
            self._merge_duplicates()

        if merge_subsets:
            LOGGER.info("Merging peptide hits that are subsets together.")
            self._merge_subsequences()

    def copy(self):
        """
        Make a copy of self.

        Returns
        -------
        :class:`DataSet<pyproteome.data_sets.DataSet>`
        """
        new = copy.copy(self)
        new.psms = new.psms.copy()
        return new

    @property
    def enrichment(self):
        """
        The str version of self.enrichments

        Returns
        -------
        str
        """
        return "+".join(self.enrichments)

    @property
    def tissue(self):
        """
        The str version of self.tissues

        Returns
        -------
        str
        """
        return "+".join(self.tissues)

    def __str__(self):
        return (
            "<pyproteome.DataSet object" +
            (
                ": " + self.name
                if self.name else
                " at " + hex(id(self))
            ) +
            ">"
        )

    def __getitem__(self, key):
        if isinstance(key, slice):
            new = self.copy()
            new.psms = new.psms[key]
            return new

        if any(
            isinstance(key, i)
            for i in [str, list, tuple, pd.Series, np.ndarray]
        ):
            return self.psms[key]

        raise TypeError(type(key))

    def _merge_duplicates(self):
        channels = list(self.channels.keys())
        agg_dict = dict((channel, sum) for channel in channels)

        agg_dict["Validated"] = all
        agg_dict["Scan Paths"] = utils.flatten_set
        agg_dict["First Scan"] = utils.flatten_set

        self.psms = self.psms.groupby(
            by=[
                "Proteins",
                "Sequence",
                "Modifications",
            ],
            sort=False,
            as_index=False,
        ).agg(agg_dict)

    def _merge_subsequences(self):
        """
        Merges petides that are a subsequence of another peptide.

        Only merges peptides that contain the same set of modifications and
        that map to the same protein(s).
        """
        psms = self.psms

        # Find all proteins that have more than one peptide mapping to them
        for index, row in psms[
            psms.duplicated(subset="Proteins", keep=False)
        ].iterrows():
            seq = row["Sequence"]

            # Then iterate over each peptide and find other non-identical
            # peptides that map to the same protein
            for o_index, o_row in psms[
                np.logical_and(
                    psms["Proteins"] == row["Proteins"],
                    psms["Sequence"] != seq,
                )
            ].iterrows():
                # If that other peptide is a subset of this peptide, rename it
                if o_row["Sequence"] in seq:
                    psms.set_value(o_index, "Sequence", seq)
                    psms.set_value(o_index, "Modifications", seq.modifications)

        # And finally group together peptides that were renamed
        self._merge_duplicates()

    def __add__(self, other):
        """
        Concatenate two data sets.

        Combines two data sets, adding together the channel values for any
        common data.

        Parameters
        ----------
        other : :class:`DataSet<pyproteome.data_sets.DataSet>`

        Returns
        -------
        :class:`DataSet<pyproteome.data_sets.DataSet>`
        """
        return merge_data([self, other])

    def normalize(self, levels, inplace=False):
        """
        Normalizes a channel to given levels.

        Divides all channel values by a given level.

        Parameters
        ----------
        levels : dict of str, float
            Mapping of channel names to normalized levels.
        inplace : bool, optional
            Modify this data set in place.

        Returns
        -------
        :class:`DataSet<pyproteome.data_sets.DataSet>`
        """
        new = self

        # Don't normalize a data set twice!
        assert not new.normalized

        if not inplace:
            new = new.copy()

        new_channels = utils.norm(self.channels)

        for key, norm_key in zip(self.channels, new_channels):
            new.psms[norm_key] = new.psms[key] / levels[key]

        new.normalized = True
        new.channels = new_channels
        new.groups = OrderedDict(
            (key, utils.norm(val))
            for key, val in self.groups.items()
        )

        new._update_snr_change()

        return new

    def dropna(self, inplace=False):
        """
        Drop any channels with NaN values.

        Parameters
        ----------
        inplace : bool, optional

        Returns
        -------
        :class:`DataSet<pyproteome.data_sets.DataSet>`
        """
        new = self

        if not inplace:
            new = new.copy()

        new.psms = new.psms.dropna(
            axis=0,
            how="any",
            subset=list(new.channels.keys()),
        )

        return new

    def filter(
        self,
        ion_score_cutoff=None, confidence_cutoff=None,
        snr_cutoff=None,
        asym_snr_cutoff=None,
        p_cutoff=None,
        fold_cutoff=None,
        asym_fold_cutoff=None,
        inplace=False,
    ):
        """
        Filters a data set.

        Parameters
        ----------
        ion_score_cutoff : int, optional
        confidence_cutoff : {"High", "Medium", "Low"}, optional
        snr_cutoff : float, optional
        asym_snr_cutoff : float, optional
        p_cutoff : float, optional
        fold_cutoff : float, optional
        asym_fold_cutoff : float, optional
        inplace : bool, optional

        Returns
        -------
        :class:`DataSet<pyproteome.data_sets.DataSet>`
        """
        if ion_score_cutoff or confidence_cutoff:
            assert self.source in ["MASCOT"]

        new = self

        if not inplace:
            new = new.copy()

        confidence_levels = []

        # MASCOT confidence levels
        if confidence_cutoff:
            confidence_levels += [confidence_cutoff]
        if confidence_cutoff in ["Medium", "Low"]:
            confidence_levels += ["High"]
        if confidence_levels in ["Low"]:
            confidence_levels += ["Medium"]

        if confidence_levels:
            new.psms = new.psms[
                new.psms["Confidence Level"].isin(confidence_levels)
            ]

        # MASCOT ion scores
        if ion_score_cutoff:
            new.psms = new.psms[
                new.psms["IonScore"] >= ion_score_cutoff
            ]

        if snr_cutoff:
            new.psms = new.psms[
                abs(new.psms["SNR"]) >= snr_cutoff
            ]

        if asym_snr_cutoff:
            if asym_snr_cutoff > 0:
                new.psms = new.psms[
                    new.psms["SNR"] >= asym_snr_cutoff
                ]
            else:
                new.psms = new.psms[
                    new.psms["SNR"] <= asym_snr_cutoff
                ]

        if p_cutoff:
            new.psms = new.psms[
                new.psms["p-value"] <= p_cutoff
            ]

        if asym_fold_cutoff:
            if asym_fold_cutoff > 1:
                new.psms = new.psms[
                    new.psms["Fold Change"] >= asym_fold_cutoff
                ]
            else:
                new.psms = new.psms[
                    new.psms["Fold Change"] <= asym_fold_cutoff
                ]

        if fold_cutoff:
            if fold_cutoff < 1:
                fold_cutoff = 1 / fold_cutoff

            new.psms = new.psms[
                np.maximum.reduce(
                    [
                        new.psms["Fold Change"],
                        1 / new.psms["Fold Change"],
                    ]
                ) >= fold_cutoff
            ]

        return new

    def _update_snr_change(self):
        """
        Update a table's SNR, Fold-Change, and p-value columns.

        Values are calculated based on changes between group_a and group_b.

        Parameters
        ----------
        psms : :class:`pandas.DataFrame`
        group_a : list of str
        group_b : list of str
        normalize : bool, optional

        Returns
        -------
        :class:`pandas.DataFrame`
        """
        groups = list(self.groups.values())

        if len(groups) == 2:
            self.psms["SNR"] = pd.Series(
                [
                    _snr(row[groups[0]], row[groups[1]])
                    for index, row in self.psms.iterrows()
                ]
            )

            self.psms["Fold Change"] = pd.Series(
                self.psms[groups[0]].mean(axis=1) /
                self.psms[groups[1]].mean(axis=1)
            )
            self.psms["p-value"] = pd.Series(
                ttest_ind(
                    self.psms[groups[0]],
                    self.psms[groups[1]],
                    axis=1,
                )[1]
            )


def merge_data(
    data_sets, name=None,
    merge_duplicates=True, merge_subsets=False,
):
    """
    Merge a list of data sets together.

    Parameters
    ----------
    data_sets : list of :class:`DataSet<pyproteome.data_sets.DataSet>`
    name : str, optional
    merge_duplicates : bool, optional
    merge_subsets : bool, optional

    Returns
    -------
    :class:`DataSet<pyproteome.data_sets.DataSet>`
    """
    assert len(data_sets) > 0

    if any(not isinstance(data, DataSet) for data in data_sets):
        raise TypeError(
            "Incompatible types: {}".format(
                [type(data) for data in data_sets]
            )
        )

    if any(data.groups != data_sets[0].groups for data in data_sets):
        raise Exception("Incongruent groups between data sets.")

    new = data_sets[0].copy()
    new.psms = pd.concat([data.psms for data in data_sets])

    if merge_duplicates:
        new._merge_duplicates()

    if merge_subsets:
        new._merge_subsequences()

    new.sets = sum(data.sets for data in data_sets)
    new.enrichments = sorted(
        set(
            enrichment
            for data in data_sets
            for enrichment in data.enrichments
        )
    )
    new.tissues = sorted(
        set(
            tissue
            for data in data_sets
            for tissue in data.tissues
        )
    )

    if new.groups:
        new._update_snr_change()

    if name:
        new.name = name

    return new


def _log_cum_fold_change(vals):
    """
    Calculate the cumulative fold change (in base-2) of values.

    Fold-change is normalized to the first element in the array.

    Parameters
    ----------
    vals : :class:`numpy.ndarray` of float

    Returns
    -------
    float
    """
    return sum(
        abs(np.log2(i / (vals[0])))
        for i in vals[1:]
        if i != 0
    )


def _snr(data1, data2):
    """
    Calculate the signal-to-noise ratio between two groups of data.

    Parameters
    ----------
    data1 : :class:`numpy.ndarray` of float
    data2 : :class:`numpy.ndarray` of float

    Returns
    -------
    float
    """
    return (data1.mean() - data2.mean()) \
        / (data1.std(ddof=1) + data2.std(ddof=1))
