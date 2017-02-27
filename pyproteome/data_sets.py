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
import numpy.ma as ma
from scipy.stats import ttest_ind

from . import loading, modification, utils


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
        dropna=False,
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

        if enrichments is None:
            enrichments = []

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

        self.intra_normalized = False
        self.inter_normalized = False
        self.sets = 1

        if dropna:
            LOGGER.info("Dropping channels with NaN values.")
            self.dropna(inplace=True)

        if filter_bad and "Confidence Level" in self.psms.columns:
            # self.psms = self.psms[
            #     self.psms["Confidence Level"].isin(["Medium", "High"])
            # ]
            self.psms = self.psms[
                self.psms["Confidence Level"].isin(["High"])
            ]

        self.update_group_changes()

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
        channels = list(self.channels.values())
        agg_dict = {}

        for channel in channels:
            if self.inter_normalized:
                self.psms[channel] = (
                    self.psms[channel] * self.psms["{}_weight".format(channel)]
                )
                agg_dict["{}_weight".format(channel)] = _nan_sum

            agg_dict[channel] = _nan_sum

        agg_dict["Validated"] = all
        agg_dict["Scan Paths"] = utils.flatten_set
        agg_dict["First Scan"] = utils.flatten_set
        agg_dict["IonScore"] = max

        self.psms = self.psms.groupby(
            by=[
                "Proteins",
                "Sequence",
                "Modifications",
            ],
            sort=False,
            as_index=False,
        ).agg(agg_dict)

        if self.inter_normalized:
            for channel in channels:
                self.psms[channel] = (
                    self.psms[channel] / self.psms["{}_weight".format(channel)]
                )

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

    def inter_normalize(self, norm_channel_name, inplace=False):
        """
        Normalize runs to one channel for inter-run comparions.

        Parameters
        ----------
        norm_channel_name : str
        inplace : bool, optional
            Modify this data set in place.

        Returns
        -------
        :class:`DataSet<pyproteome.data_sets.DataSet>`
        """

        # Don't normalize a data set twice!
        if self.inter_normalized:
            return self

        new = self
        norm_channel = self.channels[norm_channel_name]

        if not inplace:
            new = new.copy()

        for new_channel, old_channel in self.channels.items():
            new.psms[new_channel] = (
                new.psms[old_channel] / new.psms[norm_channel]
            )
            new.psms["{}_weight".format(new_channel)] = new.psms[old_channel]
            del new.psms[old_channel]

        new.inter_normalized = True
        new.channels = OrderedDict([(key, key) for key in new.channels.keys()])
        new.groups = self.groups.copy()

        new.update_group_changes()

        return new

    def normalize(self, levels, inplace=False):
        """
        Normalize channels to given levels for intra-run comparisons.

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
        assert not new.intra_normalized

        if not inplace:
            new = new.copy()

        new_channels = utils.norm(self.channels)

        for key, norm_key in zip(
            self.channels.values(),
            new_channels.values(),
        ):
            new.psms[norm_key] = new.psms[key] / levels[key]

        new.intra_normalized = True
        new.channels = new_channels
        new.groups = self.groups.copy()

        new.update_group_changes()

        return new

    def dropna(self, how="any", inplace=False):
        """
        Drop any channels with NaN values.

        Parameters
        ----------
        how : str, optional
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
            how=how,
            subset=list(new.channels.values()),
        )

        return new

    def filter(
        self,
        ion_score_cutoff=None, confidence_cutoff=None,
        p_cutoff=None,
        fold_cutoff=None,
        asym_fold_cutoff=None,
        mod_types=None,
        inplace=False,
    ):
        """
        Filters a data set.

        Parameters
        ----------
        ion_score_cutoff : int, optional
        confidence_cutoff : {"High", "Medium", "Low"}, optional
        p_cutoff : float, optional
        fold_cutoff : float, optional
        asym_fold_cutoff : float, optional
        mod_types : list of tuple of str, str, optional
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

        if mod_types:
            new.psms = modification.filter_mod_types(new.psms, mod_types)

        new.psms.reset_index(inplace=True, drop=True)

        return new

    def get_groups(self, group_a=None, group_b=None):
        """
        Get channels associated with two groups.

        Parameters
        ----------
        group_a : str or list of str, optional
        group_b : str or list of str, optional

        Returns
        -------
        groups : list of str
        labels : list of str
        """
        groups = [
            val
            for key, val in self.groups.items()
            if any(chan in self.channels for chan in val)
        ]
        labels = [
            key
            for key, val in self.groups.items()
            if any(chan in self.channels for chan in val)
        ]

        if group_a is None:
            label_a = labels[0]
            group_a = groups[0]
        elif isinstance(group_a, str):
            label_a = group_a
            group_a = self.groups[group_a]
        else:
            label_a = ", ".join(group_a)
            group_a = [
                group
                for i in group_a
                for group in self.groups[i]
            ]

        if group_b is None:
            label_b = labels[1]
            group_b = groups[1]
        elif isinstance(group_b, str):
            label_b = group_b
            group_b = self.groups[group_b]
        else:
            label_b = ", ".join(group_b)
            group_b = [
                group
                for i in group_b
                for group in self.groups[i]
            ]

        return (group_a, group_b), (label_a, label_b)

    def update_group_changes(self, group_a=None, group_b=None):
        """
        Update a table's Fold-Change, and p-value columns.

        Values are calculated based on changes between group_a and group_b.

        Parameters
        ----------
        psms : :class:`pandas.DataFrame`
        group_a : str or list, optional
        group_b : str or list, optional

        Returns
        -------
        :class:`pandas.DataFrame`
        """
        (group_a, group_b), _ = self.get_groups(
            group_a=group_a,
            group_b=group_b,
        )
        channels_a = [
            self.channels[i]
            for i in group_a
            if i in self.channels
        ]

        channels_b = [
            self.channels[i]
            for i in group_b
            if i in self.channels
        ]

        if channels_a and channels_b:
            self.psms["Fold Change"] = pd.Series(
                np.nanmean(self.psms[channels_a], axis=1) /
                np.nanmean(self.psms[channels_b], axis=1)
            )

            pvals = ttest_ind(
                self.psms[channels_a],
                self.psms[channels_b],
                axis=1,
                nan_policy="omit",
            )[1]

            if pvals is ma.masked:
                pvals = np.nan

            self.psms["p-value"] = pd.Series(pvals)
        else:
            self.psms["Fold Change"] = np.nan
            self.psms["p-value"] = np.nan


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

    if (
        any(data.channels != data_sets[0].channels for data in data_sets) and
        not all(data.inter_normalized for data in data_sets)
    ):
        raise Exception(
            "Cannot compare across runs without inter-run normalization"
        )

    new = data_sets[0].copy()
    new.psms = pd.concat([data.psms for data in data_sets])

    new.channels = OrderedDict()
    for data in data_sets:
        for key, val in data.channels.items():
            if key in new.channels:
                if new.channels[key] == val:
                    continue
                else:
                    raise Exception("Merging incompatible channels")

            new.channels[key] = val

    if merge_duplicates:
        new._merge_duplicates()

    if merge_subsets:
        new._merge_subsequences()

    new.sets = sum(data.sets for data in data_sets)
    new.enrichments = sorted(
        set(
            enrichment
            for data in data_sets
            if data.enrichments
            for enrichment in data.enrichments
        )
    )

    new.tissues = sorted(
        set(
            tissue
            for data in data_sets
            if data.tissues
            for tissue in data.tissues
        )
    )

    new.update_group_changes()

    if name:
        new.name = name

    return new


def _nan_sum(lst):
    if all(np.isnan(i) for i in lst):
        return np.nan
    else:
        return np.nansum(lst)
