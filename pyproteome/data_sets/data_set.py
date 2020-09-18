'''
This module provides functionality for manipulating proteomics data sets.

Functionality includes merging data sets and interfacing with attributes in a
structured format.
'''

# Built-ins
from __future__ import absolute_import, division

from collections import OrderedDict
import copy
import logging
import os
import re
import warnings
from itertools import chain
from functools import partial

# Core data analysis libraries
import pandas as pd
import pingouin as pg
import numpy as np
import numpy.ma as ma
from scipy.stats import ttest_ind, pearsonr, spearmanr, variation

from . import modification, protein, sequence, constand

import pyproteome as pyp


LOGGER = logging.getLogger('pyproteome.data_sets')

#:
DEFAULT_FILTER_BAD = dict(
    ion_score=15,
    isolation=30,
    median_quant=1.5e3,
    q=0.05,
)
'''
Default parameters for filtering data sets.

Selects all ions with an ion score > 15, isolation interference < 50,
median quantification signal > 1e3, and optional false-discovery q-value <
0.05.
'''

DATA_SET_COLS = [
    'Proteins',
    'Sequence',
    'Modifications',
    'Validated',
    'Confidence Level',
    'Ion Score',
    'q-value',
    'Isolation Interference',
    'Missed Cleavages',
    'Ambiguous',
    'Charges',
    'Masses',
    'RTs',
    'Intensities',
    'Raw Paths',
    'Scan Paths',
    'Scan',
    'Fold Change',
    'p-value',
]
'''
Columns available in DataSet.psms.

Note that this does not include columns for quantification or weights.
'''


class DataSet:
    '''
    Class that encompasses a proteomics data set. Data sets can be initialized
    by calling this class's constructor directly, or using
    :func:`.load_all_data`.

    Includes peptide-spectrum matches, quantification info, and mappings
    between channels, samples, and sample groups.

    Data sets are automatically loaded, filtered, and merged by default. See
    :any:`DEFAULT_FILTER_BAD<pyproteome.data_sets.data_set.DEFAULT_FILTER_BAD>`
    for default filtering parameters. See
    :func:`merge_duplicates<pyproteome.data_sets.data_set.DataSet.merge_duplicates>`
    for info on how multiple peptide-spectrum matches are integrated together.

    Attributes
    ----------
    search_name : str, optional
        Name of the search file this data set was loaded from.
    psms : :class:`pandas.DataFrame`
        Contains at least 'Proteins', 'Sequence', and 'Modifications' columns.
    channels : dict of str, str
        Maps label channel to sample name.
    groups : dict of str, list of str
        Maps groups to list of sample names. The primary group is considered
        as the first in this sequence.
    cmp_groups : list of list of str
        List of groups that are being compared.
    name : str
        Name of this data set.
    levels : dict or str, float
        Peptide levels used for normalization.
    sets : int
        Number of sets merged into this data set.
    '''
    def __init__(
        self,
        name='',
        psms=None,
        search_name=None,
        channels=None,
        groups=None,
        cmp_groups=None,
        fix_channel_names=True,
        dropna=False,
        pick_best_psm=True,
        constand_norm=False,
        merge_duplicates=True,
        filter_bad=True,
        check_raw=True,
        skip_load=False,
        skip_logging=False,
    ):
        '''
        Initializes a data set.

        Parameters
        ----------
        name : str, optional
            Name of the data set, by default this is also the .msf file that
            is used to load peptide-spectrum matches.
        psms : :class:`pandas.DataFrame`, optional
            Prebuilt dataframe of peptide-spectrum matches.
        search_name : str, optional
            Read psms from MASCOT / Discoverer data files.
        channels : dict of (str, str), optional
            Ordered dictionary mapping sample names to quantification channels
            (i.e. {'X': '126', 'Y': '127', 'Z': '128', 'W': '129'})
        groups : dict of str, list of str, optional
            Ordered dictionary mapping sample names to larger groups
            (i.e. {'WT': ['X', 'Y'], 'Diseased': ['W', 'Z']})
        cmp_groups : list of list of str
            cmp_groups that are passed to :func:`.DataSet.norm_cmp_groups()`.
        fix_channel_names : bool, optional
        dropna : bool, optional
            Drop scans that have any channels with missing quantification
            values.
        pick_best_psm : bool, optional
            Select the peptide sequence for each scan that has the highest
            MASCOT ion score. (i.e. ['pSTY': 5, 'SpTY': 10, 'STpY': 20] =>
            'STpY')
        constand_norm : bool, optional
        merge_duplicates : bool, optional
            Merge scans that have the same peptide sequence into one peptide,
            summing the quantification channel intensities to give a weighted
            estimate of relative abundances.
        filter_bad : bool or dict, optional
            Remove peptides that do not have a 'High' confidence score from
            ProteomeDiscoverer.
        check_raw : bool, optional
        skip_load : bool, optional
            Just initialize the structure, don't load any data.
        skip_logging : bool, optional
            Don't log any information.
        '''
        if search_name is None:
            search_name = name

        if search_name and os.path.splitext(search_name)[1] == '':
            search_name += '.msf'

        self.channels = channels.copy() if channels else OrderedDict()
        self.groups = groups.copy() if groups else OrderedDict()
        self.cmp_groups = cmp_groups or None
        self.group_a, self.group_b = None, None

        species, lst = set(), []

        if search_name and not skip_load:
            self.psms, species, lst = pyp.loading.load_psms(
                search_name,
                pick_best_psm=pick_best_psm,
            )
            for col in DATA_SET_COLS:
                assert col in self.psms.columns
        else:
            self.psms = pd.DataFrame(
                columns=DATA_SET_COLS + list(self.channels.values()),
            )

        if psms:
            self.add_peptide(psms)

        self.name = name
        self.levels = None
        self.search_name = search_name
        self.species = species

        self.intra_normalized = False
        self.inter_normalized = False
        self.sets = 1

        if fix_channel_names:
            self.fix_channel_names()

        if check_raw:
            self.check_raw()

        if dropna:
            LOGGER.info(
                '{}: Dropping channels with NaN values.'.format(self.name)
            )
            self.dropna(inplace=True)

        if filter_bad is True:
            filter_bad = DEFAULT_FILTER_BAD.copy()

            if pd.isnull(self.psms['q-value']).all() and 'q' in filter_bad:
                del filter_bad['q']

            if pd.isnull(self.psms[list(self.channels.values())]).all().all():
                del filter_bad['median_quant']

        if not skip_logging:
            self.log_stats()

        if filter_bad:
            LOGGER.info(
                '{}: Filtering peptides using: {}'
                .format(self.name, filter_bad)
            )
            self.filter(
                filter_bad,
                inplace=True,
            )
            self.log_stats()

        if pick_best_psm and (
            not search_name or
            os.path.splitext(search_name)[1] != '.msf' or
            any(i for i in lst)
        ):
            LOGGER.info(
                '{}: Picking peptides with best ion score for each scan.'
                .format(self.name)
            )
            self._pick_best_psm()

        if constand_norm:
            channels = list(self.channels.values())

            # Display quant distribution
            pyp.levels.get_channel_levels(self)

            # Save channel intensities as weights for data integration
            for channel in channels:
                weight = '{}_weight'.format(channel)
                self.psms[weight] = self.psms[channel]

            # Run CONSTANd normalization
            constand.constand(self, inplace=True)

        if merge_duplicates:
            LOGGER.info(
                '{}: Merging duplicate peptide hits together.'
                .format(self.name)
            )
            self.merge_duplicates(inplace=True)

        if cmp_groups:
            self.norm_cmp_groups(cmp_groups, inplace=True)

        self.update_group_changes()

        if not skip_logging:
            self.log_stats()

    def copy(self):
        '''
        Make a copy of self.

        Returns
        -------
        ds : :class:`.DataSet`
        '''
        new = copy.copy(self)

        new.psms = new.psms.copy()
        new.channels = new.channels.copy()
        new.groups = new.groups.copy()
        new.species = new.species.copy()

        return new

    def fix_channel_names(self):
        '''
        Correct quantification channel names to those present in the search
        file.

        i.e. from 130_C to 130C (or vice versa).
        '''
        for key, val in list(self.channels.items()):
            if '_' in val:
                new_val = val.replace('_', '')
            else:
                new_val = val[:-1] + '_' + val[-1]

            if (
                val not in self.psms.columns and
                new_val in self.psms.columns
            ):
                LOGGER.info(
                    '{}: Correcting {}: {} to {}'.format(
                        self.name,
                        key, val, new_val,
                    )
                )
                self.channels[key] = new_val

    @property
    def samples(self):
        '''
        Get a list of sample names in this data set.

        Returns
        -------
        list of str
        '''
        return self.get_samples()

    def get_samples(self, groups=None):
        '''
        Get a list of sample names in this data set.

        Parameters
        ----------
        groups : optional, list of (list of str)

        Returns
        -------
        list of str
        '''
        if groups is None:
            groups = self.cmp_groups or [list(self.groups.keys())]
        
        channel_names = [
            sample_name
            for lst in groups
            for group in lst
            for sample_name in self.groups[group]
            if sample_name in self.channels
        ]
        # Remove duplicates
        channel_names = list(OrderedDict.fromkeys(channel_names))
        return channel_names

    def __str__(self):
        return (
            '<pyproteome.DataSet object' +
            (
                ': ' + self.name
                if self.name else
                ' at ' + hex(id(self))
            ) +
            '>'
        )

    def __getitem__(self, key):
        if isinstance(key, slice):
            new = self.copy()
            new.psms = new.psms[key]
            return new

        if any(
            isinstance(key, i)
            for i in [str, list, set, tuple, pd.Series, np.ndarray]
        ):
            return self.psms[key]

        raise TypeError(type(key))

    def __setitem__(self, key, val):
        if any(
            isinstance(key, i)
            for i in [str, list, set, tuple, pd.Series, np.ndarray]
        ):
            self.psms[key] = val
        else:
            raise TypeError(type(key))

    @property
    def shape(self):
        '''
        Get the size of a data set in (rows, columns) format.

        Returns
        -------
        shape : tuple of (int, int)
        '''
        return self.psms.shape

    def _pick_best_psm(self):
        reject_mask = np.zeros(self.shape[0], dtype=bool)

        for index, row in self.psms.iterrows():
            hits = np.logical_and(
                self.psms['Scan'] == row['Scan'],
                self.psms['Sequence'] != row['Sequence'],
            )

            if 'Rank' in self.psms.columns:
                better = self.psms['Rank'] < row['Rank']
            else:
                better = self.psms['Ion Score'] > row['Ion Score']

            hits = np.logical_and(hits, better)

            if hits.any():
                reject_mask[index] = True

        self.psms = self.psms[~reject_mask].reset_index(drop=True)

    def merge_duplicates(self, inplace=False):
        '''
        Merge together all duplicate peptides. New quantification values are
        calculated from a weighted sum of each channel's values.

        Parameters
        ----------
        inplace : bool, optional
            Modify the data set in place, otherwise create a copy and return
            the new object.

        Returns
        -------
        ds : :class:`.DataSet`
        '''
        new = self

        if not inplace:
            new = new.copy()

        if new.shape[0] < 1:
            return new

        channels = list(new.channels.values())
        agg_dict = OrderedDict()

        for channel in channels:
            weight = '{}_weight'.format(channel)
            std = '{}_std'.format(channel)
            mean = '{}_mean'.format(channel)
            cv = '{}_cv'.format(channel)

            if weight in new.psms.columns:
                new.psms[channel] *= new.psms[weight]
                agg_dict[weight] = 'sum'

            agg_dict[channel] = 'sum'

            if cv not in new.psms.columns:
                new.psms[std] = new.psms[channel]
                new.psms[mean] = new.psms[channel]

                agg_dict[std] = 'std'
                agg_dict[mean] = 'mean'
            else:
                agg_dict[cv] = 'mean'

        agg_dict['Modifications'] = _first
        agg_dict['Missed Cleavages'] = _first
        agg_dict['Validated'] = all

        agg_dict['Scan Paths'] = pyp.utils.flatten_set
        agg_dict['Raw Paths'] = pyp.utils.flatten_set

        agg_dict['Ambiguous'] = all

        agg_dict['Masses'] = pyp.utils.flatten_set
        agg_dict['Charges'] = pyp.utils.flatten_set
        agg_dict['Intensities'] = pyp.utils.flatten_set
        agg_dict['RTs'] = pyp.utils.flatten_set

        agg_dict['Scan'] = pyp.utils.flatten_set
        agg_dict['Ion Score'] = max
        agg_dict['q-value'] = min
        agg_dict['Confidence Level'] = partial(
            max,
            key=lambda x: ['Low', 'Medium', 'High'].index(x),
        )
        agg_dict['Isolation Interference'] = min

        new.psms = new.psms.groupby(
            by=[
                'Proteins',
                'Sequence',
            ],
            sort=False,
            as_index=False,
        ).agg(agg_dict).reset_index(drop=True)

        for channel in channels:
            weight = '{}_weight'.format(channel)
            std = '{}_std'.format(channel)
            mean = '{}_mean'.format(channel)
            cv = '{}_cv'.format(channel)

            if cv not in new.psms.columns:
                new.psms[cv] = new.psms[std] / new.psms[mean]
                del new.psms[std]
                del new.psms[mean]

            if weight in new.psms.columns:
                new.psms[channel] = (
                    new.psms[channel] / new.psms[weight]
                )

            new.psms[channel] = new.psms[channel].replace(0, np.nan)

        return new

    def merge_subsequences(self, inplace=False):
        '''
        Merges petides that are a subsequence of another peptide.
        (i.e. SVYTEIKIHK + SVYTEIK -> SVYTEIK)

        Only merges peptides that contain the same set of modifications and
        that map to the same protein(s).

        Parameters
        ----------
        inplace : bool, optional

        Returns
        -------
        ds : :class:`.DataSet`
        '''
        new = self

        if not inplace:
            new = new.copy()

        new['__sort__'] = new['Sequence'].apply(len)
        new.psms = new.psms.sort_values('__sort__', ascending=True)

        # Find all proteins that have more than one peptide mapping to them
        for index, row in new.psms[
            new.psms.duplicated(subset='Proteins', keep=False)
        ].iterrows():
            seq = row['Sequence']

            # Then iterate over each peptide and find other non-identical
            # peptides that map to the same protein
            for o_index, o_row in new.psms[
                np.logical_and(
                    new.psms['Proteins'] == row['Proteins'],
                    new.psms['Sequence'] != seq,
                )
            ].iterrows():
                # If that other peptide is a subset of this peptide, rename it
                if o_row['Sequence'] in seq:
                    cols = [
                        'Sequence',
                        'Modifications',
                        'Missed Cleavages',
                    ]
                    new.psms.at[o_index, cols] = row[cols]


        del new.psms['__sort__']

        # And finally group together peptides that were renamed
        return new.merge_duplicates(inplace=inplace)

    def __add__(self, other):
        '''
        Concatenate two data sets.

        Combines two data sets, adding together the channel values for any
        common data.

        Parameters
        ----------
        other : :class:`.DataSet`

        Returns
        -------
        ds : :class:`.DataSet`
        '''
        return merge_data([self, other])

    def rename_channels(self, inplace=False):
        '''
        Rename all columns names for quantification channels to sample names.
        (i.e. '126' => 'Mouse #1').

        Parameters
        ----------
        inplace : bool, optional

        Returns
        -------
        ds : :class:`.DataSet`
        '''
        new = self

        if not inplace:
            new = new.copy()

        for new_channel, old_channel in new.channels.items():
            if new_channel != old_channel:
                new_weight = '{}_weight'.format(new_channel)

                if (
                    new_channel in new.psms.columns or
                    new_weight in new.psms.columns
                ):
                    raise Exception(
                        'Channel {} already exists, cannot rename to it'
                        .format(new_channel)
                    )

                new.psms[new_channel] = new.psms[old_channel]
                del new.psms[old_channel]

                old_weight = '{}_weight'.format(old_channel)

                if old_weight in new.psms.columns:
                    new.psms[new_weight] = new.psms[old_weight]
                    del new.psms[old_weight]

        new.channels = OrderedDict([
            (key, key) for key in new.channels.keys()
        ])

        return new

    def inter_normalize(
        self, 
        norm_channels=None, 
        other=None, 
        inplace=False,
    ):
        '''
        Normalize runs to one channel for inter-run comparions.

        Parameters
        ----------
        other : :class:`.DataSet`, optional
            Second data set to normalize quantification values against,
            using a common normalization channels.
        norm_channels : list of str, optional
            Normalization channels to use for cross-run normalization.
        inplace : bool, optional
            Modify this data set in place.

        Returns
        -------
        ds : :class:`.DataSet`
        '''
        assert (
            norm_channels is not None or
            other is not None
        )

        new = self

        if not inplace:
            new = new.copy()

        new = new.rename_channels(inplace=inplace)

        if norm_channels is None:
            norm_channels = set(new.channels).intersection(other.channels)

        # If there are no common channels, do not change anything
        if len(norm_channels) == 0:
            return new

        # Filter norm channels to include only those in other data set
        norm_channels = [
            chan
            for chan in norm_channels
            if not other or chan in other.channels
        ]
        for channel in new.channels.values():
            weight = '{}_weight'.format(channel)

            if weight not in new.psms.columns:
                new.psms[weight] = (
                    new.psms[channel] *
                    (100 - new.psms['Isolation Interference']) / 100
                )

        # Calculate the mean normalization signal from each shared channel
        new_mean = new.psms[norm_channels].mean(axis=1)

        # Drop values for which there is no normalization data
        new.psms = new.psms[~new_mean.isnull()].reset_index(drop=True)

        if other:
            merge = pd.merge(
                new.psms,
                other.psms,
                on=[
                    # 'ProteinStr',
                    # 'PeptideStr',
                    'Proteins',
                    'Sequence',
                    'Modifications',
                ],
                how='left',
                suffixes=('_self', '_other'),
            )
            assert merge.shape[0] == new.shape[0]

            self_mean = merge[
                ['{}_self'.format(i) for i in norm_channels]
            ].mean(axis=1)

            other_mean = merge[
                ['{}_other'.format(i) for i in norm_channels]
            ].mean(axis=1)

            for channel in new.channels.values():
                vals = merge[
                    channel
                    if channel in merge.columns else
                    '{}_self'.format(channel)
                ]

                # Set scaling factor to 1 where other_mean is None
                cp = other_mean.copy()
                cp[other_mean.isnull()] = (
                    self_mean[other_mean.isnull()]
                )

                if self_mean.any():
                    vals *= cp / self_mean

                assert new.shape[0] == vals.shape[0]
                new.psms[channel] = vals

        new.update_group_changes()
        new.inter_normalized = True

        return new

    def normalize(self, lvls, inplace=False):
        '''
        Normalize channels to given levels for intra-run comparisons.

        Divides all channel values by a given level.

        Parameters
        ----------
        lvls : dict of str, float
            Mapping of channel names to normalized levels. All quantification
            values for each channel are divided by the normalization factor.
        inplace : bool, optional
            Modify this data set in place.

        Returns
        -------
        ds : :class:`.DataSet`
        '''
        new = self

        # Don't normalize a data set twice!
        assert not new.intra_normalized

        if not inplace:
            new = new.copy()

        if hasattr(lvls, 'levels'):
            if not lvls.levels:
                lvls.levels = pyp.levels.get_channel_levels(lvls)[1]

            lvls = lvls.levels

        new_channels = pyp.utils.norm(self.channels)

        for key, norm_key in zip(
            self.channels.values(),
            new_channels.values(),
        ):
            new.psms[norm_key] = new.psms[key] / lvls[key]
            del new.psms[key]

        new.intra_normalized = True
        new.channels = new_channels
        new.groups = self.groups.copy()

        new.update_group_changes()

        return new

    def check_raw(self):
        '''
        Checks that all raw files referenced in search data can be found in
        :const:`pyproteome.paths.MS_RAW_DIR`.

        Returns
        -------
        found_all : bool
        '''
        try:
            raw_dir = [
                i.lower()
                for i in os.listdir(pyp.paths.MS_RAW_DIR)
            ]
        except OSError:
            raw_dir = []

        found_all = True

        for raw in sorted(
            pyp.utils.flatten_set(
                row['Raw Paths']
                for _, row in self.psms.iterrows()
            )
        ):
            if raw.lower() not in raw_dir:
                LOGGER.warning(
                    '{}: Unable to locate raw file for {}'
                    .format(self.name, raw)
                )
                found_all = False

        return found_all

    def dropna(
        self,
        columns=None,
        how=None,
        thresh=None,
        groups=None,
        inplace=False,
    ):
        '''
        Drop any channels with NaN values.

        Parameters
        ----------
        columns : list of str, optional
        how : str, optional
        groups : list of str, optional
            Only drop rows with NaN in columns within groups.
        inplace : bool, optional

        Returns
        -------
        ds : :class:`.DataSet`
        '''
        new = self

        if not inplace:
            new = new.copy()

        if columns is None:
            columns = list(new.channels.values())

        if how is None and thresh is None:
            how = 'any'

        if groups is not None:
            columns = [
                new.channels[col]
                for group in groups
                for col in new.groups[group]
                if col in new.channels
            ]

        new.psms = new.psms.dropna(
            axis=0,
            how=how,
            thresh=thresh,
            subset=columns,
        ).reset_index(drop=True)

        return new

    def add_peptide(self, inserts):
        '''
        Manually add a peptide or list of peptides to a data set.

        Parameters
        ----------
        insert : dict or list of dict

        Examples
        --------
        
        This example demostrates how to manually insert a peptide that was manually
        validated by the user::

            prots = data_sets.protein.Proteins(
                proteins=(
                    data_sets.protein.Protein(
                        accession='Q920G3',
                        gene='Siglec5',
                        description='Sialic acid-binding Ig-like lectin 5',
                        full_sequence=(
                            'MRWAWLLPLLWAGCLATDGYSLSVTGSVTVQEGLCVFVACQVQYPNSKGPVFGYWFREGA'
                            'NIFSGSPVATNDPQRSVLKEAQGRFYLMGKENSHNCSLDIRDAQKIDTGTYFFRLDGSVK'
                            'YSFQKSMLSVLVIALTEVPNIQVTSTLVSGNSTKLLCSVPWACEQGTPPIFSWMSSALTS'
                            'LGHRTTLSSELNLTPRPQDNGTNLTCQVNLPGTGVTVERTQQLSVIYAPQKMTIRVSWGD'
                            'DTGTKVLQSGASLQIQEGESLSLVCMADSNPPAVLSWERPTQKPFQLSTPAELQLPRAEL'
                            'EDQGKYICQAQNSQGAQTASVSLSIRSLLQLLGPSCSFEGQGLHCSCSSRAWPAPSLRWR'
                            'LGEGVLEGNSSNGSFTVKSSSAGQWANSSLILSMEFSSNHRLSCEAWSDNRVQRATILLV'
                            'SGPKVSQAGKSETSRGTVLGAIWGAGLMALLAVCLCLIFFTVKVLRKKSALKVAATKGNH'
                            'LAKNPASTINSASITSSNIALGYPIQGHLNEPGSQTQKEQPPLATVPDTQKDEPELHYAS'
                            'LSFQGPMPPKPQNTEAMKSVYTEIKIHKC'
                        ),
                    ),
                ),
            )
            seq = data_sets.sequence.extract_sequence(prots, 'SVyTEIK')

            mods = data_sets.modification.Modifications(
                mods=[
                    data_sets.modification.Modification(
                        rel_pos=0,
                        mod_type='TMT10plex',
                        nterm=True,
                        sequence=seq,
                    ),
                    data_sets.modification.Modification(
                        rel_pos=2,
                        mod_type='Phospho',
                        sequence=seq,
                    ),
                    data_sets.modification.Modification(
                        rel_pos=6,
                        mod_type='TMT10plex',
                        sequence=seq,
                    ),
                ],
            )

            seq.modifications = mods

            ckh_sigf1_py_insert = {
                'Proteins': prots,
                'Sequence': seq,
                'Modifications': mods,
                '126':  1.46e4,
                '127N': 2.18e4,
                '127C': 1.88e4,
                '128N': 4.66e3,
                '128C': 6.70e3,
                '129N': 7.88e3,
                '129C': 1.03e4,
                '130N': 7.28e3,
                '130C': 2.98e3,
                '131':  6.01e3,
                'Validated': True,
                'First Scan': {23074},
                'Raw Paths': {'2019-04-24-CKp25-SiglecF-1-py-SpinCol-col189.raw'},
                'Scan Paths': {'CK-7wk-H1-pY'},
                'IonScore': 30,
                'Isolation Interference': 0,
            }
            ds.add_peptide([ckh_sigf1_py_insert])
        '''
        defaults = {
            'Proteins': protein.Proteins(),
            'Sequence': sequence.Sequence(),
            'Modifications': modification.Modifications(),
            'Validated': False,
            'Confidence Level': 'High',
            'Ion Score': 100,
            'q-value': 0,
            'Isolation Interference': 0,
            'Missed Cleavages': 0,
            'Ambiguous': False,
            'Charges': set(),
            'Masses': set(),
            'RTs': set(),
            'Intensities': set(),
            'Raw Paths': set(),
            'Scan Paths': set(),
            'Scan': set(),
            'Fold Change': np.nan,
            'p-value': np.nan,
        }

        if isinstance(inserts, dict):
            inserts = [inserts]

        for insert in inserts:
            for key, val in defaults.items():
                if key not in insert:
                    insert[key] = val

            for col in DATA_SET_COLS:
                assert col in insert

            self.psms = self.psms.append(pd.Series(insert), ignore_index=True)

    def filter(
        self,
        filters=None,
        inplace=False,
        **kwargs
    ):
        '''
        Filters a data set.

        Parameters
        ----------
        filters : list of dict or dict, optional
            List of filters to apply to data set. Filters are also pulled from
            `kwargs` (see below).
        inplace : bool, optional
            Perform the filter on self, other create a copy and return the
            new object.

        Notes
        -----
        These parameters filter your data set to only include peptides that
        match a given attribute. For example::

            >>> data.filter(mod='Y', p=0.01, fold=2)

        This function interprets both the argument filter and python `kwargs`
        magic. The three functions are all equivalent::

            >>> data.filter(p=0.01)
            >>> data.filter([{'p': 0.01}])
            >>> data.filter({'p': 0.01})

        Filter parameters can be one of any below:

        ================    ===================================================
        Name                Description
        ================    ===================================================
        series              Use a pandas series (data.psms[series]).
        fn                  Use data.psms.apply(fn).
        group_a             Calculate p / fold change values from group_a.
        group_b             Calculate p / fold change values from group_b.
        ambiguous           Include peptides with ambiguous PTMs if true,
                            filter them out if false.
        confidence          Discoverer's peptide confidence (High|Medium|Low).
        ion_score           MASCOT's ion score.
        isolation           Discoverer's isolation inference.
        missed_cleavage     Missed cleaves <= cutoff.
        median_quant        Median quantification signal >= cutoff.
        median_cv           Median coefficient of variation <= cutoff.
        p                   p-value < cutoff.
        q                   q-value < cutoff.
        asym_fold           Change > val if cutoff > 1 else Change < val.
        fold                Change > cutoff or Change < 1 / cutoff.
        motif               Filter for motif.
        protein             Filter for protein or list of proteins.
        protein             Filter for protein or list of proteins.
        accession           Filter for protein or list of UniProt accessions.
        sequence            Filter for sequence or list of sequences.
        mod                 Filter for modifications.
        only_validated      Use rows validated by CAMV.
        any                 Use rows that many any filter.
        inverse             Use all rows that are rejected by a filter.
        rename              Change the new data sets name to a new value.
        ================    ===================================================

        Returns
        -------
        ds : :class:`.DataSet`
        '''
        new = self

        if filters is None:
            filters = []

        if filters and not isinstance(filters, (list, tuple)):
            filters = [filters]

        if kwargs:
            filters += [kwargs]

        if not inplace:
            new = new.copy()

        filters = copy.deepcopy(filters)

        confidence = {
            'High': ['High'],
            'Medium': ['Medium', 'High'],
            'Low': ['Low', 'Medium', 'High'],
        }

        fns = {
            'series': lambda val, psms:
            val,

            'fn': lambda val, psms:
            psms.apply(val, axis=1),

            'ambiguous': lambda val, psms:
            psms['Ambiguous'] == val,

            'confidence': lambda val, psms:
            psms['Confidence Level'].isin(confidence[val]),

            'ion_score': lambda val, psms:
            psms['Ion Score'] >= val,

            'isolation': lambda val, psms:
            psms['Isolation Interference'] <= val,

            'missed_cleavage': lambda val, psms:
            psms['Missed Cleavages'] <= val,

            'median_quant': lambda val, psms:
            np.nan_to_num(
                np.nanmedian(
                    psms[
                        [chan for chan in new.channels.values()]
                    ],
                    axis=1,
                )
            ) >= val,

            'median_cv': lambda val, psms:
            np.nan_to_num(
                np.nanmedian(
                    psms[
                        [
                            '{}_CV'.format(chan)
                            for chan in new.channels.values()
                        ]
                    ],
                    axis=1,
                )
            ) <= val,

            'p': lambda val, psms:
            (~psms['p-value'].isnull()) &
            (psms['p-value'] <= val),

            'q': lambda val, psms:
            (~psms['q-value'].isnull()) &
            (psms['q-value'] <= val),

            'asym_fold': lambda val, psms:
            (~psms['Fold Change'].isnull()) &
            (
                psms['Fold Change'] >= val
                if val > 1 else
                psms['Fold Change'] <= val
            ),

            'fold': lambda val, psms:
            (~psms['Fold Change'].isnull()) &
            (
                psms['Fold Change'].apply(
                    lambda x:
                    x if x > 1 else (1 / x if x else x)
                ) >= (val if val > 1 else 1 / val)
            ),

            'motif': lambda val, psms:
            psms['Sequence'].apply(
                lambda x:
                any(
                    (
                        pyp.motifs.motif.Motif(val)
                        if isinstance(val, str) else
                        val
                    ).match(nmer)
                    for nmer in pyp.motifs.generate_n_mers(
                        x,
                        mods=f.get('mod', None),
                    )
                )
            ),

            'protein': lambda val, psms:
            psms['Proteins'].apply(
                lambda x: bool(set(val).intersection(x.genes))
            )
            if isinstance(val, (list, set, tuple, pd.Series)) else
            psms['Proteins'] == val.strip(),

            'accession': lambda val, psms:
            psms['Proteins'].apply(
                lambda x: bool(set(val).intersection(x.accessions))
            )
            if isinstance(val, (list, set, tuple, pd.Series)) else
            psms['Proteins'].apply(
                lambda x: any([i == val.strip() for i in x.accessions])
            ),

            'sequence': lambda val, psms:
            psms['Sequence'].apply(lambda x: any(i in x for i in val))
            if isinstance(val, (list, set, tuple, pd.Series)) else
            psms['Sequence'] == val.strip(),

            'mod': lambda val, psms:
            psms['Modifications'].apply(
                lambda x: bool(list(x.get_mods(val).skip_labels()))
            ),

            'only_validated': lambda val, psms:

            psms['Validated'] == val,

            'scan_paths': lambda val, psms:
            psms['Scan Paths']
            .apply(lambda x: any(i in val for i in x))
        }

        with warnings.catch_warnings():
            warnings.filterwarnings(
                'ignore',
                r'All-NaN (slice|axis) encountered',
            )

            for f in filters:
                group_a = f.pop('group_a', None)
                group_b = f.pop('group_b', None)

                if group_a or group_b:
                    new.update_group_changes(
                        group_a=group_a,
                        group_b=group_b,
                    )

                inverse = f.pop('inverse', False)
                rename = f.pop('rename', None)

                if rename is not None:
                    new.name = rename

                if f.pop('any', False):
                    if new.shape[0] < 1:
                        return new

                    mask = pd.Series(
                        [inverse] * new.shape[0],
                        index=new.psms.index,
                    )

                    for key, val in f.items():
                        # Skip filtering if psms is empty
                        f_mask = fns[key](val, new.psms)

                        if inverse:
                            f_mask = ~f_mask

                        mask |= f_mask

                    assert mask.shape[0] == new.shape[0]
                    new.psms = new.psms.loc[mask].reset_index(drop=True)
                else:
                    for key, val in f.items():
                        # Skip filtering if psms is empty
                        if new.shape[0] < 1:
                            continue

                        mask = fns[key](val, new.psms)

                        if inverse:
                            mask = ~mask

                        assert mask.shape[0] == new.shape[0]

                        new.psms = new.psms.loc[mask].reset_index(drop=True)

        new.psms.reset_index(inplace=True, drop=True)

        return new

    def get_groups(self, group_a=None, group_b=None):
        '''
        Get channels associated with two groups.

        Parameters
        ----------
        group_a : str or list of str, optional
        group_b : str or list of str, optional

        Returns
        -------
        samples : list of str
        labels : list of str
        groups : tuple of (str or list of str)
        '''
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

        group_a = group_a or self.group_a
        group_b = group_b or self.group_b

        if group_a is None:
            label_a = labels[0] if labels else None
            samples_a = groups[0] if groups else []
        elif isinstance(group_a, str):
            label_a = group_a
            samples_a = self.groups[group_a]
        else:
            group_a = [
                group
                for group in group_a
                if any(
                    sample in self.channels
                    for sample in self.groups[group]
                )
            ]
            label_a = ', '.join(group_a)
            samples_a = [
                sample
                for i in group_a
                for sample in self.groups[i]
                if sample in self.channels
            ]

        if group_b is None:
            label_b = labels[1] if labels[1:] else None
            samples_b = groups[1] if groups[1:] else []
        elif isinstance(group_b, str):
            label_b = group_b
            samples_b = self.groups[group_b]
        else:
            group_b = [
                group
                for group in group_b
                if any(
                    sample in self.channels
                    for sample in self.groups[group]
                )
            ]
            label_b = ', '.join(group_b)
            samples_b = [
                sample
                for i in group_b
                for sample in self.groups[i]
                if sample in self.channels
            ]

        return (samples_a, samples_b), (label_a, label_b), (group_a, group_b)

    def update_group_changes(self, group_a=None, group_b=None):
        '''
        Update a DataSet's Fold Change, and p-value for each peptide using the give two-group comparison.

        Values are calculated based on changes between group_a and group_b. p-values are calculated as a 2-sample t-test.

        Parameters
        ----------
        psms : :class:`pandas.DataFrame`
        group_a : str or list of str, optional
            Single or multiple groups to use for fold change numerator.
        group_b : str or list of str, optional
            Single or multiple groups to use for fold change denominator.
        '''
        (group_a, group_b), _, (self.group_a, self.group_b) = self.get_groups(
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

        if channels_a and channels_b and self.shape[0] > 0:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', category=RuntimeWarning)
                self.psms['Fold Change'] = pd.Series(
                    np.nanmean(self.psms[channels_a], axis=1) /
                    np.nanmean(self.psms[channels_b], axis=1),
                    index=self.psms.index,
                )

                pvals = ttest_ind(
                    self.psms[channels_a],
                    self.psms[channels_b],
                    axis=1,
                    nan_policy='omit',
                )[1]

                if ma.is_masked(pvals):
                    pvals = ma.fix_invalid(pvals, fill_value=np.nan)
                elif isinstance(pvals, float):
                    pvals = [pvals]
                elif pvals.shape == ():
                    pvals = [pvals]

                self.psms['p-value'] = pd.Series(pvals, index=self.psms.index)
        else:
            self.psms['Fold Change'] = np.nan
            self.psms['p-value'] = np.nan

    def norm_cmp_groups(self, cmp_groups, ctrl_groups=None, inplace=False):
        '''
        Normalize between groups in a list. This can be used to compare data
        sets that have comparable control groups.

        Channnels within each list of groups are normalized to the mean of the
        group's channels.

        Parameters
        ----------
        cmp_groups : list of list of str
            List of groups to be normalized to each other.
            i.e. [('CK Hip', 'CK-p25 Hip'), ('CK Cortex', 'CK-p25 Cortex')]
        ctrl_groups : list of str, optional
            List of groups to use for baseline normalization. If not set,
            the first group from each comparison will be used.
        inplace : bool, optional
            Modify the data set in place, otherwise create a copy and return
            the new object.

        Examples
        --------
        >>> channels = ['a', 'b', 'c', 'd']
        >>> groups = {i: [i] for i in channels}
        >>> ds = data_sets.DataSet(channels=channels, groups=groups)
        >>> ds.add_peptide({'a': 1000, 'b': 500, 'c': 100, 'd': 25})
        >>> ds = ds.norm_cmp_groups([['a', 'b'], ['c', 'd']])
        >>> ds.data
        {'a': 1, 'b': 0.5, 'c': 1, 'd': .25}

        Returns
        -------
        ds : :class:`.DataSet`
        '''
        new = self

        if not inplace:
            new = new.copy()

        if new.cmp_groups:
            raise Exception('cmp_groups normalization already set')

        new.cmp_groups = [
            tuple(i) for i in cmp_groups
        ]

        # if len(cmp_groups) < 2:
        #     return new

        # assert not any(
        #     group in [
        #         i
        #         for o_groups in cmp_groups[ind + 1:]
        #         for i in o_groups
        #     ]
        #     for ind, groups in enumerate(cmp_groups)
        #     for group in groups
        # )
        for groups in cmp_groups:
            channels = [
                [
                    new.channels[name]
                    for name in new.groups[group]
                    if name in new.channels
                ]
                for group in groups
            ]
            group_names = [
                group
                for ind, group in enumerate(groups)
            ]

            if not any([len(i) > 0 for i in channels]):
                continue

            vals = new[
                [chan for chan_group in channels for chan in chan_group]
            ]

            if ctrl_groups:
                norm_vals = vals[[
                    chan
                    for i in ctrl_groups
                    if i in group_names
                    for chan in channels[group_names.index(i)]
                ]].median(axis=1)
            else:
                norm_vals = vals[
                    channels[0]
                ].median(axis=1)

            vals = vals.apply(lambda col: col / norm_vals)

            new.psms[
                [chan for chan_group in channels for chan in chan_group]
            ] = vals

        return new

    def log_stats(self):
        '''
        Log statistics information about peptides contained in the data
        set. This information includes total numbers, phospho-specificity,
        modification ambiguity, completeness of labeling, and missed
        cleavage counts.
        '''
        data_p = self.filter(mod='Phospho')
        data_no_p = self.filter(mod='Phospho', inverse=True)
        data_pst = self.filter(mod=[('S', 'Phospho'), ('T', 'Phospho')])
        data_py = self.filter(mod=[('Y', 'Phospho')])

        LOGGER.info('{}: Data Set Statistics:'.format(self.name))

        LOGGER.info(
            (
                '{}: -- {} pY - {} pST - {} total ({:.0%} phospho specificity)'
            ).format(
                self.name,
                len(data_py.psms),
                len(data_pst.psms),
                len(self.psms),
                len(data_p.psms) / max([len(self.psms), 1]),
            )
        )
        LOGGER.info(
            (
                '{}: -- {} unique proteins'
            ).format(
                self.name,
                len(self.genes),
            )
        )

        if len(data_p.psms) > 0:
            per_prot_quant = (
                len(set(data_p.genes).intersection(data_no_p.genes)) /
                len(data_p.psms)
            )

            LOGGER.info(
                (
                    '{}: -- {:.0%} of phosphopeptides have '
                    'protein quantification'
                ).format(self.name, per_prot_quant)
            )

        per_amb = (
            data_p.filter(ambiguous=True).shape[0] /
            max([data_p.shape[0], 1])
        )

        LOGGER.info(
            (
                '{}: -- {:.0%} of phosphopeptides have an ambiguous assignment'
            ).format(self.name, per_amb)
        )

        labeled = self.filter(
            fn=lambda x:
            x['Sequence'].is_labeled
        ).shape[0]

        underlabeled = self.filter(
            fn=lambda x:
            x['Sequence'].is_underlabeled
        ).shape[0]

        per_lab = labeled / max([self.shape[0], 1])
        per_under_lab = underlabeled / max([self.shape[0], 1])

        LOGGER.info(
            (
                '{}: -- {:.0%} labeled - {:.0%} underlabeled'
            ).format(self.name, per_lab, per_under_lab)
        )

        mmc = self.psms['Missed Cleavages'].mean()

        LOGGER.info(
            (
                '{}: -- {:.1f} mean missed cleavages'
            ).format(self.name, mmc)
        )

    @property
    def genes(self):
        '''
        Get all uniprot gene names occuring in this data set.

        Returns
        -------
        list of str
        '''
        return sorted(
            set(
                gene
                for i in self.psms['Proteins']
                for gene in i.genes
            )
        )

    @property
    def phosphosites(self):
        '''
        Get a list of all unique phosphosites identified in a data set.

        Returns
        -------
        list of str

        Examples
        --------
            >>> ds.phosphosites
            ['Siglec5 pY561', 'Stat3 pY705', 'Inpp5d pY868']
        '''
        return sorted(
            set(
                self.psms.apply(
                    lambda x:
                    '{0}{1}{2}'.format(
                        pyp.utils.get_name(x['Proteins']),
                        ' : '
                        if len(x['Modifications'].get_mods('Phospho')) > 0 else
                        '',
                        re.sub(
                            r'(\d+)',
                            r'\1',
                            x['Modifications'].get_mods('Phospho').__str__(
                                prot_index=0,
                                show_mod_type=False,
                            ),
                        ),
                    ),
                    axis=1,
                ).values
            )
        )

    @property
    def accessions(self):
        '''
        Get all uniprot accessions occuring in this data set.

        Returns
        -------
        list of str

        Examples
        --------
            >>> ds.accessions
            ['P42227', 'Q920G3', 'Q9ES52']
        '''
        return sorted(
            set(
                gene
                for i in self.psms['Proteins']
                for gene in i.accessions
            )
        )

    @property
    def data(self):
        '''
        Get the quantification data for all samples and peptides in a data set.

        Returns
        -------
        df : :class:`pandas.DataFrame`
        '''
        return self.get_data()

    def get_data(self, groups=None, mods=None, short_name=False):
        '''
        Get the quantification data for all samples and peptides in a data set, with
        more customizable options.

        Parameters
        ----------
        groups : list of str, optional
            Select samples from a given list of groups, otherwise select all samples.
        mods : str or list of str, optional
            Option passed to :func:`.modification.Modification.__str__`.
        short_name : bool, optional
            Use the short abbreviation of a gene name (using :func:`pyproteome.utils.get_name`).
            Otherwise use the long version.

        Returns
        -------
        df : :class:`pandas.DataFrame`
        '''
        if groups is None:
            channel_names = self.samples

        d = self.psms[
            [
                self.channels[chan]
                for chan in channel_names
            ]
        ]
        d.index = self.psms.apply(
            lambda row:
            '{} {} - {}'.format(
                (
                    pyp.utils.get_name(row['Proteins'])
                    if short_name else
                    ' / '.join(row['Proteins'].genes)[:16]
                ),
                (
                    row['Modifications'].get_mods(mods)
                    if mods else
                    row['Modifications']
                ).__str__(prot_index=0),
                row['Sequence'],
            ),
            axis=1,
        )

        return d

    @property
    def intensity_data(self, groups=None, norm_cmp=False):
        '''
        Get the quantification data for all samples and peptides in a data set.

        Parameters
        ----------
        norm_cmp : bool, optional

        Returns
        -------
        df : :class:`pandas.DataFrame`
        '''
        if groups is None:
            channel_names = self.samples

        d = self.psms[
            [
                self.channels[chan] + '_weight'
                for chan in channel_names
            ]
        ]
        d.index = self.psms.apply(
            lambda row:
            '{} {} - {}'.format(
                ' / '.join(row['Proteins'].genes)[:16],
                row['Modifications'].__str__(prot_index=0),
                row['Sequence'],
            ),
            axis=1,
        )
        d.columns = [i.rsplit('_', 1)[0] for i in d.columns]

        return d


def load_all_data(
    chan_mapping=None,
    group_mapping=None,
    loaded_fn=None,
    norm_mapping=None,
    merge_mapping=None,
    merged_fn=None,
    kw_mapping=None,
    merge_only=True,
    replace_norm=True,
    **kwargs
):
    '''
    Load, normalize, and merge all data sets found in
    :const:`pyproteome.paths.MS_SEARCHED_DIR`.

    Parameters
    ----------
    chan_mapping : dict, optional
    group_mapping : dict, optional
    loaded_fn : func, optional
    norm_mapping : dict, optional
    merge_mapping : dict, optional
    merged_fn : func, optional
    kw_mapping : dict of (str, dict)
    merge_only : bool, optional
    replace_norm : bool, optional
        If true, only keep the normalized version of a data set. Otherwise
        return both normalized and unnormalized version.
    kwargs : dict
        Any extra arguments are passed directly to :class:`.DataSet` during
        initialization.

    Returns
    -------
    datas : dict of str, :class:`.DataSet`

    Examples
    --------
    This example demostrates how to automatically load, filter, normalize, and
    together several data sets::

        ckh_channels = OrderedDict(
            [
                ('3130 CK Hip',     '126'),
                ('3131 CK-p25 Hip', '127'),
                ('3145 CK-p25 Hip', '128'),
                ('3146 CK-p25 Hip', '129'),
                ('3148 CK Hip',     '130'),
                ('3157 CK Hip',     '131'),
            ]
        )
        ckx_channels = OrderedDict(
            [
                ('3130 CK Cortex',     '126'),
                ('3131 CK-p25 Cortex', '127'),
                ('3145 CK-p25 Cortex', '128'),
                ('3146 CK-p25 Cortex', '129'),
                ('3148 CK Cortex',     '130'),
                ('3157 CK Cortex',     '131'),
            ]
        )
        ckp25_groups = OrderedDict(
            [
                (
                    'CK',
                    [
                        '3130 CK Hip',
                        '3148 CK Hip',
                        '3157 CK Hip',
                        '3130 CK Cortex',
                        '3148 CK Cortex',
                        '3157 CK Cortex',
                    ],
                ),
                (
                    'CK-p25',
                    [
                        '3131 CK-p25 Hip',
                        '3145 CK-p25 Hip',
                        '3146 CK-p25 Hip',
                        '3131 CK-p25 Cortex',
                        '3145 CK-p25 Cortex',
                        '3146 CK-p25 Cortex',
                    ],
                ),
            ]
        )
        # With search data located as follows:
        #   Searched/
        #       CK-H1-pY.msf
        #       CK-H1-pST.msf
        #       CK-H1-Global.msf
        #       CK-X1-pY.msf
        #       CK-X1-pST.msf
        #       CK-X1-Global.msf
        # Load each data set, normalized to its respective global proteome analysis:
        datas = data_sets.load_all_data(
            chan_mapping={
                'CK-H': ckh_channels,
                'CK-X': ckx_channels,
            },
            # Normalize pY, pST, and Global runs to each sample's global data
            norm_mapping={
                'CK-H1': 'CK-H1-Global',
                'CK-X1': 'CK-X1-Global',
            ]),
            # Merge together normalized hippocampus and cortex runs
            merge_mapping={
                'CK Hip': ['CK-H1-pY', 'CK-H1-pST', 'CK-H1-Global'],
                'CK Cortex': ['CK-X1-pY', 'CK-X1-pST', 'CK-X1-Global'],
                'CK All': ['CK Hip', 'CK Cortex'],
            },
            groups=ckp25_groups,
        )

        # Alternatively, load each data set, using CONSTANd normalization:
        data_sets.constand.DEFAULT_CONSTAND_COL = 'kde'
        datas = data_sets.load_all_data(
            chan_mapping={
                'CK-H': ckh_channels,
                'CK-X': ckx_channels,
            },
            norm_mapping='constand',
            # Merge together normalized hippocampus and cortex runs
            merge_mapping={
                'CK Hip': ['CK-H1-pY', 'CK-H1-pST', 'CK-H1-Global'],
                'CK Cortex': ['CK-X1-pY', 'CK-X1-pST', 'CK-X1-Global'],
                'CK All': ['CK Hip', 'CK Cortex'],
            },
            groups=ckp25_groups,
        )
    '''
    chan_mapping = chan_mapping or {}
    group_mapping = group_mapping or {}
    norm_mapping = norm_mapping or {}
    merge_mapping = merge_mapping or {}
    kw_mapping = kw_mapping or {}

    datas = OrderedDict()

    for f in sorted(os.listdir(pyp.paths.MS_SEARCHED_DIR)):
        name, ext = os.path.splitext(f)

        if ext not in ['.msf']:
            continue

        if (
            merge_mapping and
            merge_only and
            name not in pyp.utils.flatten_set(list(merge_mapping.values()))
        ):
            continue

        kws = kwargs.copy()
        kws.update(kw_mapping.get(name, {}))

        # Run CONSTANd normalization before data filtering and integration
        if 'constand_norm' not in kws and norm_mapping == 'constand':
            kws['constand_norm'] = True

        chan = kws.pop('channels', None)
        group = kws.pop('groups', None)

        for key, val in chan_mapping.items():
            if name.startswith(key):
                chan = val
                break

        for key, val in group_mapping.items():
            if name.startswith(key):
                group = val
                break

        datas[name] = DataSet(
            name=name,
            channels=chan,
            groups=group,
            **kws
        )

        if loaded_fn:
            datas[name] = loaded_fn(name, datas[name])

    datas, mapped_names = norm_all_data(
        datas,
        norm_mapping,
        replace_norm=replace_norm,
    )
    datas = merge_all_data(
        datas, merge_mapping,
        mapped_names=mapped_names,
        merged_fn=merged_fn,
    )

    return datas


def norm_all_data(
    datas,
    norm_mapping,
    replace_norm=True,
    inplace=True,
):
    '''
    Normalize all data sets.

    Parameters
    ----------
    datas : dict of (str, :class:`.DataSet`)
    norm_mapping : dict of (str, str)
    replace_norm : bool, optional
    inplace : bool, optional
        Modify datas object inplace.

    Returns
    -------
    datas : dict of (str, :class:`.DataSet`)
    mapped_names : dict of (str, str)
    '''
    mapped_names = OrderedDict()

    if not inplace:
        datas = datas.copy()

    constand_norm = norm_mapping in ['constand']

    if norm_mapping in ['self', 'constand']:
        norm_mapping = {
            name: name + ('' if replace_norm else '-norm')
            for name in datas.keys()
            if not name.endswith('-norm')
        }

    for name, data in list(datas.items()):
        if name.endswith('-norm'):
            continue

        for key, val in norm_mapping.items():
            if not name.startswith(key):
                continue

            mapped_names[name] = '{}-norm'.format(name)

            if not constand_norm:
                data = data.normalize(datas[val], inplace=replace_norm)
            else:
                # data = constand.constand(
                #     data, 
                #     inplace=replace_norm,
                # )
                pass

            data.name += '-norm'

            datas[mapped_names[name]] = data

            if replace_norm:
                del datas[name]

            break

    return datas, mapped_names


def merge_all_data(
    datas,
    merge_mapping,
    mapped_names=None,
    merged_fn=None,
    inplace=True,
):
    '''
    Merge together multiple data sets.

    Parameters
    ----------
    datas : dict of (str, :class:`.DataSet`)
    merge_mapping : dict of (str, list of str)
    mapped_names : dict of (str, str), optional
    merged_fn : func, optional
    inplace : bool, optional
        Modify datas object inplace.

    Returns
    -------
    datas : dict of (str, :class:`.DataSet`)
    '''
    if not inplace:
        datas = datas.copy()

    if mapped_names is None:
        mapped_names = {}

    rmap = {val: key for key, val in mapped_names.items()}

    if merge_mapping:
        for name in datas.keys():
            if not any(
                name in vals or
                rmap.get(name, name) in vals or
                name == key
                for key, vals in merge_mapping.items()
            ):
                LOGGER.warning('Unmerged data: {}'.format(name))

    for key, vals in merge_mapping.items():
        if not any([mapped_names.get(val, val) in datas for val in vals]):
            continue

        datas[key] = merge_data(
            [
                datas[mapped_names.get(val, val)]
                for val in vals
                if mapped_names.get(val, val) in datas
            ],
            name=key,
        )

        if merged_fn:
            datas[key] = merged_fn(key, datas[key])

    return datas


def merge_data(
    data_sets,
    name=None,
    norm_channels=None,
    merge_duplicates=True,
):
    '''
    Merge a list of data sets together.

    Parameters
    ----------
    data_sets : list of :class:`.DataSet`
    name : str, optional
    norm_channels : dict of (str, str)
    merge_duplicates : bool, optional

    Returns
    -------
    ds : :class:`.DataSet`
    '''
    new = DataSet(
        name=name,
        skip_logging=True,
        skip_load=True,
        filter_bad=False,
        constand_norm=False,
    )

    if len(data_sets) < 1:
        return new

    # if any(not isinstance(data, DataSet) for data in data_sets):
    #     raise TypeError(
    #         'Incompatible types: {}'.format(
    #             [type(data) for data in data_sets]
    #         )
    #     )

    for index, data in enumerate(data_sets):
        data = data.rename_channels()

        # Update new.groupss
        for group, samples in data.groups.items():
            if group not in new.groups:
                new.groups[group] = samples
                continue

            new.groups[group] += [
                sample
                for sample in samples
                if sample not in new.groups[group]
            ]

        if not data.inter_normalized:
            # Normalize data sets to their common channels
            inter_norm = norm_channels

            if not inter_norm:
                inter_norm = set(data.channels)

                if index > 0:
                    inter_norm = inter_norm.intersection(new.channels)
                elif len(data_sets) > 1:
                    inter_norm = inter_norm.intersection(
                        data_sets[1].channels
                    )

            data = data.inter_normalize(
                other=new if index > 0 else None,
                norm_channels=inter_norm,
            )

        for key, val in data.channels.items():
            assert new.channels.get(key, val) == val

            if key not in new.channels:
                new.channels[key] = val

        new.psms = _concat(
            [new.psms, data.psms],
        ).reset_index(drop=True)

        if merge_duplicates:
            new.merge_duplicates(inplace=True)

    new.sets = sum(data.sets for data in data_sets)
    new.inter_normalized = True

    new.group_a = next(
        (i.group_a for i in data_sets if i.group_a is not None),
        None,
    )
    new.group_b = next(
        (i.group_b for i in data_sets if i.group_b is not None),
        None,
    )
    new.cmp_groups = sorted(
        set(
            grp
            for i in data_sets
            if i.cmp_groups is not None
            for grp in i.cmp_groups
        )
    ) or None
    new.species = set(
        species
        for i in data_sets
        for species in i.species
    )

    new.update_group_changes()

    if name:
        new.name = name

    new.log_stats()

    return new

def _unfirst(x):
    return x.values[0]

def _first(x):
    if not all(i == x.values[0] for i in x.values):
        LOGGER.warning(
            'Mismatch between peptide data: \'{}\' not in {}'
            .format(
                x.values[0],
                [str(i) for i in x.values[1:]],
            )
        )

    return x.values[0]

def merge_proteins(ds, inplace=False, fn=None):
    '''
    Merge together all peptides mapped to the same protein. Maintains the
    first available peptide and calculates the median quantification value
    for each protein across all of its peptides.

    Parameters
    ----------
    ds : :class:`.DataSet`
    inplace : bool, optional

    Returns
    -------
    ds : :class:`.DataSet`
    '''
    new = ds

    if not inplace:
        new = new.copy()

    if len(new.psms) < 1:
        return new

    if fn is None:
        fn = 'median'

    channels = list(new.channels.values())
    agg_dict = {}

    for channel in channels:
        weight = '{}_weight'.format(channel)

        if weight in new.psms.columns:
            # XXX: Use weight corresponding to median channel value?
            # (not median weight)
            agg_dict[weight] = fn

        agg_dict[channel] = fn


    agg_dict['Sequence'] = _unfirst
    agg_dict['Modifications'] = _unfirst
    agg_dict['Missed Cleavages'] = _unfirst
    agg_dict['Validated'] = all

    agg_dict['Scan Paths'] = _unfirst
    agg_dict['Raw Paths'] = _unfirst

    agg_dict['Ambiguous'] = _unfirst

    agg_dict['Masses'] = _unfirst
    agg_dict['Charges'] = _unfirst
    agg_dict['Intensities'] = _unfirst
    agg_dict['RTs'] = _unfirst

    agg_dict['Scan'] = _unfirst
    agg_dict['Ion Score'] = max
    agg_dict['q-value'] = min
    agg_dict['Confidence Level'] = partial(
        max,
        key=lambda x: ['Low', 'Medium', 'High'].index(x),
    )
    agg_dict['Isolation Interference'] = min
    agg_dict['Unique Peptides'] = len

    new.psms['Unique Peptides'] = np.nan

    new.psms = new.psms.groupby(
        by=[
            'Proteins',
        ],
        sort=False,
        as_index=False,
    ).agg(agg_dict)

    new.update_group_changes()

    return new


def _cv(x):
    x = variation(x, nan_policy='omit')
    return np.nan if ma.is_masked(x) else x


def update_correlation(ds, corr, metric='spearman', min_periods=5):
    '''
    Update a table's Fold-Change, and p-value columns.

    Values are calculated based on changes between group_a and group_b.

    Parameters
    ----------
    ds : :class:`.DataSet`
    corr : :class:`pandas.Series`
    metric : str, optional
    min_periods : int, optional

    Returns
    -------
    ds : :class:`.DataSet`
    '''
    ds = ds.copy()

    corr = pd.to_numeric(
        corr[(~corr.isnull()) & (corr.index.isin(ds.psms.columns))]
    )

    metric = {
        'spearman': partial(spearmanr, b=corr, nan_policy='omit'),
        'pearson':
        lambda x: 
        pearsonr(
            x[~x.isnull()], 
            y=corr[~x.isnull()],
        ),
    }[metric]

    # Avoids errors with correlation of 2 points
    cols = list(corr.index)
    ds = ds.dropna(thresh=3, columns=cols)

    df = ds.psms[cols].apply(pd.to_numeric)

    vals = df.apply(metric, axis=1, result_type='expand')

    vals[
        (~df.isnull()).sum(axis=1) < min_periods
    ] = np.nan

    ds.psms[['Fold Change', 'p-value']] = vals

    return ds


def _pair_corr(row, x_val, y_val):
    pcorr_df = pd.DataFrame([
        pd.to_numeric(row[x_val.index.tolist()]),
        x_val,
        y_val
    ]).T
    pcorr_df.columns = ['row', 'x', 'y']

    x_pair_corr = pcorr_df.pairwise_corr(
        covar='y',
        method='spearman',
    )
    y_pair_corr = pcorr_df.pairwise_corr(
        covar='x',
        method='spearman',
    )
    x_corr = x_pair_corr['r'].iloc[0]
    y_corr = y_pair_corr['r'].iloc[0]

    x_p = x_pair_corr['p-unc'].iloc[0]
    y_p = y_pair_corr['p-unc'].iloc[0]
    
    return pd.Series(
        (x_corr, y_corr, x_p, y_p),
        index=[
        'Correlation-x',
        'Correlation-y',
        'p-value-x',
        'p-value-y',
    ])


def update_pairwise_corr(ds, x_val, y_val, inplace=False):
    if not inplace:
        ds = ds.copy()
    
    vals = ds.psms.apply(
        partial(_pair_corr, x_val=x_val, y_val=y_val),
        axis=1,
        result_type='expand',
    )
    ds[vals.columns.tolist()] = vals
    
    return ds


def _concat(dfs):
    cols = []

    for frame in dfs:
        for col in frame.columns:
            if col not in cols:
                cols.append(col)

    df_dict = dict.fromkeys(cols, [])

    def _fast_flatten(input_list):
        return list(chain.from_iterable(input_list))

    for col in cols:
        # Flatten and save to df_dict
        df_dict[col] = list(
            chain.from_iterable(
                frame[col]
                if col in frame.columns else
                [np.nan] * frame.shape[0]
                for frame in dfs
            )
        )

    df = pd.DataFrame.from_dict(df_dict)[cols]

    return df
