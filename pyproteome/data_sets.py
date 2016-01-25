
# Built-ins
import copy
import logging

# Core data analysis libraries
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind

from . import loading, math, utils


LOGGER = logging.getLogger("pyproteome.DataSet")


class DataSet:
    """
    Class that encompasses a proteomics data set.

    Includes peptide list, scan list, channels, groups, quantitative
    phenotypes, etc.

    Attributes
    ----------
    psms : pandas.DataFrame
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
    enrichment : str
    tissue : str
    sets : int
        Number of sets merged into this data set.
    """
    def __init__(
        self, channels,
        psms=None,
        mascot_name=None,
        camv_name=None,
        groups=None, phenotypes=None,
        name="", enrichment="", tissue="",
        dropna=True,
        camv_slices=1,
    ):
        """
        Initializes a data set.

        Parameters
        ----------
        channels : dict of str, str
        psms : pandas.DataFrame, optional
        mascot_name : str, optional
        camv_name : str, optional
        groups : dict of str, list of str, optional
        phenotypes : dict of str, (dict of str, float), optional
        name : str, optional
        enrichment : str, optional
        tissue : str, optional
        dropna : bool, optional
        camv_slices : int, optional
        """
        assert psms is not None or \
            mascot_name is not None or \
            camv_name is not None

        self.source = "unknown"

        if mascot_name:
            psms = loading.load_mascot_psms(
                mascot_name,
                camv_slices=1,
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
        self.enrichment = enrichment
        self.tissue = tissue

        self.normalized = False
        self.sets = 1

        if dropna:
            LOGGER.info("Dropping channels with NaN values.")
            self.dropna(inplace=True)

    def copy(self):
        new = copy.copy(self)
        new.psms = new.psms.copy()
        return new

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

        if isinstance(key, str):
            return self.psms[key]

        raise TypeError

    def _merge_psms(self):
        channels = list(self.channels.keys())

        if self.normalized:
            channels += utils.norm(channels)

        agg_dict = dict((channel, sum) for channel in channels)
        self.psms = self.psms.groupby(
            [
                "Proteins",
                "Sequence",
                "Modifications",
            ],
            sort=False,
            as_index=False,
        ).agg(agg_dict)

    def __add__(self, other):
        """
        Concatenate two data sets.

        Combines two data sets, adding together the channel values for any
        common data.

        Parameters
        ----------
        other : pyproteome.DataSet

        Returns
        -------
        pyproteome.DataSet
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
        pyproteome.DataSet
        """
        new = self

        # Don't normalize a data set twice!
        assert not new.normalized

        if not inplace:
            new = new.copy()

        for key, norm_key in zip(self.channels, utils.norm(self.channels)):
            new.psms[norm_key] = new.psms[key] / levels[key]

        new.normalized = True

        new._merge_psms()
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
        pyproteome.DataSet
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
        pyproteome.DataSet
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
                    new.psms["SNR"] >= snr_cutoff
                ]
            else:
                new.psms = new.psms[
                    new.psms["SNR"] <= snr_cutoff
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
        psms : pandas.DataFrame
        group_a : list of str
        group_b : list of str
        normalize : bool, optional

        Returns
        -------
        pandas.DataFrame
        """
        groups = list(self.groups.values())

        if self.normalized:
            groups = utils.norm(groups)

        if len(groups) == 2:
            self.psms["SNR"] = pd.Series(
                [
                    math.snr(row[groups[0]], row[groups[1]])
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


def merge_data(data_sets):
    """
    Merge a list of data sets together.

    Parameters
    ----------
    data_sets : list of pyproteome.DataSet

    Returns
    -------
    pyproteome.DataSet
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
    new._merge_psms()
    new.sets = sum(data.sets for data in data_sets)

    if new.groups:
        new._update_snr_change()

    return new
