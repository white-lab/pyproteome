'''
Plot calculated levels of a given sequence across channels or groups.
'''

from __future__ import division

# Built-ins
from collections import OrderedDict
import itertools
import logging

# Core data analysis libraries
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import ttest_ind

import pyproteome as pyp

LOGGER = logging.getLogger('pyproteome.plot')


def plot(
    data,
    title=None,
    ax=None,
    log_2=True,
    box=True,
):
    '''
    Plot the levels of a sequence across multiple channels.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.DataSet`
    title : str, optional
    figsize : tuple of (int, int), optional

    Returns
    -------
    figs : list of :class:`matplotlib.figure.Figure`
    '''

    cmp_groups = data.cmp_groups or [list(data.groups.keys())]

    channel_names = [
        channel_name
        for groups in cmp_groups
        for group in groups
        for channel_name in data.groups[group]
        if channel_name in data.channels
    ]
    channel_names = list(OrderedDict.fromkeys(channel_names))
    channels = [
        data.channels[channel_name]
        for channel_name in channel_names
    ]

    figures = []

    for _, row in data.psms.iterrows():
        seq = str(row['Sequence'])

        values = row[channels]
        mask = ~pd.isnull(row[channels])
        names = pd.Series(channel_names, index=values.index)[mask]
        values = values[mask]

        values = values / (
            values[[
                data.channels[i]
                for grp in cmp_groups
                for i in data.groups[grp[0]]
                if i in data.channels and data.channels[i] in values.index
            ]].median()
            if len(cmp_groups) > 0 else
            values[0]
        )

        if ax is None:
            fig, ax_i = plt.subplots(
                figsize=(len(values) / 2, 6 / 2),
            )
        else:
            ax_i = ax

        df = pd.DataFrame(
            [
                (
                    name,
                    val,
                    [
                        group
                        for groups in cmp_groups
                        for group in groups
                        if name in data.groups[group]
                    ][0],
                )
                for name, val in zip(names, values)
            ],
            columns=('name', 'val', 'group'),
        )

        if log_2:
            df['val'] = df['val'].apply(np.log2)

        (sns.barplot if box else sns.barplot)(
            x='name',
            y='val',
            hue='group',
            data=df,
            ax=ax_i,
            dodge=False,
        )
        ax_i.set_xticklabels(
            ax_i.get_xticklabels(),
            rotation=45,
            horizontalalignment='right',
        )
        ax_i.set_xlabel('')
        ax_i.get_legend().set_title('')

        ax_i.axhline(
            np.log2(1) if log_2 else 1,
            linestyle=':',
            color='#4C4D4F',
            # alpha=.5,
        )

        mod_str = row['Modifications'].__str__(prot_index=0)

        ax_i.set_title(
            title
            if title else
            '{} ({}{})'.format(
                seq,
                pyp.utils.get_name(row['Proteins'])[:20],
                (' ' + mod_str) if mod_str else '',
            ),
        )

        ylabel = 'Intensity'
        # if data.cmp_groups:
        # else:
        #     ylabel = (
        #         'Cummulative Intensity' +
        #         (' (Normalized)' if data.intra_normalized else '')
        #     )

        ax_i.set_ylabel(
            ylabel,
        )

        figures.append((ax_i.get_figure(), ax_i))

    return figures


def _gen_groups(cmp_groups):
    ret = []

    for i in cmp_groups:
        lst = list(itertools.combinations(i, 2))
        for ind, i in enumerate(lst, start=1):
            ret.append(i + (ind != len(lst),))

    return ret


def plot_group(
    data,
    cmp_groups=None,
    cmp_groups_star=None,
    title=None,
    ax=None,
    box=True,
    show_p=True,
    show_ns=False,
    log_2=True,
    offset_frac=20,
    title_mods=None,
    size=4,
    y_max=None,
    p_ha='center',
    cmap='cool',
):
    '''
    Plot the levels of a sequence across each group.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.DataSet`
    cmp_groups : list of tuple, optional
    cmp_groups_star : list of tuple, optional
    title : str, optional
    ax : :class:`matplotlib.axes.Axes`, optional
    box : bool, optional
    show_p : bool, optional
    show_ns : bool, optional
    log_2 : bool, optional
    offset_frac : float, optional
    title_mods : list of str, optional
    size : float, optional
    y_max : float, optional
    p_ha : str, optional
    cmap : str, optional

    Returns
    -------
    figs : list of :class:`matplotlib.figure.Figure`
    '''
    if cmp_groups is None:
        cmp_groups = data.cmp_groups or [list(data.groups.keys())]

    figures = []

    for _, row in data.psms.iterrows():
        values = []

        for groups in cmp_groups:
            groups = [i for i in groups if i in data.groups]
            group_vals = pd.Series([
                row[[
                    data.channels[name]
                    for name in data.groups[group]
                    if name in data.channels
                ]]
                for group in groups
            ], index=groups, dtype=object)

            group_vals = pd.Series([
                group[(~pd.isnull(group)) & (group > 0)]
                for group in group_vals
            ], index=group_vals.index, dtype=object)

            group_vals = group_vals[
                group_vals.apply(lambda x: x.shape[0] > 0)
            ]

            # Check normalization group is not null and at least one other
            # group of channels is not null
            if (
                group_vals.shape[0] < 1
                # or (
                #     len(cmp_groups) > 1
                #     # and groups[0] not in group_vals.index
                # )
                # or all([
                #     group not in group_vals.index
                #     for group in groups[1:]
                # ])
            ):
                continue

            normalize = group_vals.iloc[0].median()

            group_vals = pd.Series([
                group / normalize
                for group in group_vals
            ], index=group_vals.index, dtype=object)

            values.append(group_vals)

        labels = [
            name
            for group in values
            for name in group.index
        ]

        if ax is None:
            _, plot_ax = plt.subplots(
                figsize=(len(labels) * .75, 4),
            )
        else:
            plot_ax = ax

        x = [
            ind
            for ind, l in enumerate(
                j
                for i in values
                for j in i.values
            )
            for k in l
        ]
        y = np.concatenate([
            j.astype(float)
            for i in values
            for j in i.values
        ])

        if log_2:
            y = np.log2(y)

        def _get_color(label):
            lst = [
                i
                for i in cmp_groups
                for ind, val in enumerate(i)
                if val == label
            ][0]
            return sns.color_palette(
                cmap, len(lst),
            ).as_hex()[lst.index(label)]

        df = pd.DataFrame(
            [
                (
                    k,
                    np.log2(k),
                    label,
                )
                for i in values
                for label, j in i.iteritems()
                for k in j.values
            ],
            columns=('y', 'log2_y', 'label'),
        )

        if box:
            kwargs = {'showfliers': False}
        else:
            kwargs = {}
        
        (sns.boxplot if box else sns.barplot)(
            x='label',
            y='log2_y' if log_2 else 'y',
            hue='label',
            palette={label: _get_color(label) for label in df['label']},
            data=df,
            ax=plot_ax,
            dodge=False,
            linewidth=.75 if box else 0,
            **kwargs
        )
        sns.swarmplot(
            x=x,
            y=y,
            color='.25',
            ax=plot_ax,
            size=size,
        )
        plot_ax.axhline(
            np.log2(1) if log_2 else 1,
            linestyle='--',
            color='#4C4D4F',
        )

        mod_str = row['Modifications'].get_mods(title_mods).__str__(prot_index=0)

        plot_ax.set_title(
            title
            if title else
            '{}{}({}{})'.format(
                row['Sequence'],
                ' ' if len(plot_ax.get_xticklabels()) > 2 else '\n',
                pyp.utils.get_name(row['Proteins'])[:20],
                (' ' + mod_str) if mod_str else '',
            ),
        )
        plot_ax.xaxis.grid(False)

        if show_p:
            if y_max is None:
                y_max_cp = y.max()
            else:
                y_max_cp = y_max

            v = [
                vals
                for group_vals in values
                for vals in group_vals
            ]

            cmp_star = cmp_groups_star
            
            if cmp_star is None:
                cmp_star = _gen_groups(cmp_groups)

            offset = y_max_cp / offset_frac / 2

            for x_ind, x in enumerate(cmp_star):
                move_offset = True

                if len(x) == 4:
                    group_a, group_b, move_offset, p_ha = x
                elif len(x) == 3:
                    group_a, group_b, move_offset = x
                else:
                    group_a, group_b = x

                if group_a not in labels or group_b not in labels:
                    continue

                index_a = labels.index(group_a)
                index_b = labels.index(group_b)

                values_a = v[index_a]
                values_b = v[index_b]

                if values_a.shape[0] < 2 or values_b.shape[0] < 2:
                    continue

                pval = ttest_ind(
                    values_a.values,
                    values_b.values,
                ).pvalue
                txt = pyp.utils.stars(pval)

                if txt != 'ns' or show_ns:
                    plot_ax.annotate(
                        '',
                        xy=(
                            index_a,
                            y_max_cp + offset,
                        ),
                        xytext=(
                            index_b,
                            y_max_cp + offset,
                        ),
                        xycoords='data',
                        textcoords='data',
                        arrowprops=dict(arrowstyle='-', ec='#000000'),
                    )
                    p_x = {
                        'left': index_a,
                        'center': np.mean([index_a, index_b]),
                        'right': index_b,
                    }.get(p_ha, index_b)
                    plot_ax.text(
                        x=p_x,
                        y=y_max_cp + offset + y_max_cp / offset_frac / 4 * (2 if txt == 'ns' else .75),
                        s=txt,
                        ha='center',
                        va='center',
                    )

                    if move_offset and x_ind != len(cmp_star) - 1:
                        offset += y_max_cp / offset_frac
        else:
            y_max_cp, offset = 0, 0


        plot_ax.set_ylim(
            bottom=plot_ax.get_ylim()[0],
            top=max([
                plot_ax.get_ylim()[1],
                y_max_cp + offset + y_max_cp / offset_frac,
            ]),
        )

        plot_ax.yaxis.set_ticks_position('left')
        plot_ax.xaxis.set_ticks_position('bottom')

        plot_ax.set_xlabel('')
        plot_ax.set_ylabel(
            '{} Signal'.format(
                'Relative' if cmp_groups else 'Cumulative',
            ),
        )
        if plot_ax.get_legend():
            plot_ax.get_legend().set_visible(False)

        if log_2:
            plot_ax.set_yticks([
                i
                for i in plot_ax.get_yticks()
                if i >= plot_ax.get_ylim()[0] and i <= plot_ax.get_ylim()[1]
            ])
            plot_ax.set_yticklabels(
                [
                    '{:.2f}'.format(np.power(2, i))
                    for i in plot_ax.get_yticks()
                ],
            )

        plot_ax.set_xticklabels(
            labels,
            rotation=45,
            horizontalalignment='right',
        )

        figures.append((plot_ax.get_figure(), plot_ax))

    return figures


def plot_together(
    data,
    cmp_groups=None,
    title=None,
    ax=None,
    show_p=True,
    log_2=True,
    cmap='cool',
):
    '''
    Plot the levels of a sequence across each group in one shared plot.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.DataSet`
    cmp_groups : list of tuple, optional
    title : str, optional
    ax : :class:`matplotlib.axes.Axes`, optional
    show_p : bool, optional
    log_2 : bool, optional
    cmap : str, optional

    Returns
    -------
    figs : list of :class:`matplotlib.figure.Figure`
    '''
    if cmp_groups is None:
        cmp_groups = data.cmp_groups or [list(data.groups.keys())]

    figures = []

    for _, row in data.psms.iterrows():
        values = []

        for groups in cmp_groups:
            groups = [i for i in groups if i in data.groups]
            group_vals = pd.Series([
                row[[
                    data.channels[name]
                    for name in data.groups[group]
                    if name in data.channels
                ]]
                for group in groups
            ], index=groups, dtype=object)

            group_vals = pd.Series([
                group[~pd.isnull(group)]
                for group in group_vals
            ], index=group_vals.index, dtype=object)

            group_vals = group_vals[
                group_vals.apply(lambda x: x.shape[0] > 0)
            ]

            # Check normalization group is not null and at least one other
            # group of channels is not null
            if (
                group_vals.shape[0] < 1 or (
                    len(cmp_groups) > 1 and
                    groups[0] not in group_vals.index
                ) or all(
                    group not in group_vals.index
                    for group in groups[1:]
                )
            ):
                continue

            normalize = group_vals.iloc[0].median()

            group_vals = pd.Series([
                group / normalize
                for group in group_vals
            ], index=group_vals.index, dtype=object)

            values.append(group_vals)

        labels = [
            name
            for group in values
            for name in group.index
        ]

        if ax is None:
            _, plot_ax = plt.subplots(
                figsize=(len(labels) * .75, 4),
            )
        else:
            plot_ax = ax

        x = [
            ind
            for ind, l in enumerate(
                j
                for i in values
                for j in i.values
            )
            for k in l
        ]
        y = np.concatenate([
            j.astype(float)
            for i in values
            for j in i.values
        ])

        if log_2:
            y = np.log2(y)

        def _get_color(label):
            lst = [
                i
                for i in cmp_groups
                for ind, val in enumerate(i)
                if val == label
            ][0]
            return sns.color_palette(
                cmap, len(lst),
            ).as_hex()[lst.index(label)]

        df = pd.DataFrame(
            [
                (
                    k,
                    np.log2(k),
                    label,
                    _get_color(label),
                )
                for i in values
                for label, j in i.iteritems()
                for k in j.values
            ],
            columns=('y', 'log2_y', 'label', 'color'),
        )
        sns.boxplot(
            x='label',
            y='log2_y' if log_2 else 'y',
            hue='color',
            data=df,
            ax=plot_ax,
            dodge=False,
            boxprops=dict(alpha=.3),
        )
        sns.swarmplot(
            x=x,
            y=y,
            color='.25',
            ax=plot_ax,
            # size=10,
        )
        plot_ax.axhline(
            np.log2(1) if log_2 else 1,
            linestyle='--',
            color='#4C4D4F',
        )

        if title:
            plot_ax.set_title(title)

        plot_ax.xaxis.grid(False)

        if show_p:
            y_max = y.max()

            v = [
                vals
                for group_vals in values
                for vals in group_vals
            ]

            for grp_set in cmp_groups:
                offset = y_max / 10

                for label_a, label_b in itertools.combinations(grp_set, 2):
                    if label_a not in labels or label_b not in labels:
                        continue

                    index_a = labels.index(label_a)
                    index_b = labels.index(label_b)

                    values_a = v[index_a]
                    values_b = v[index_b]

                    if values_a.shape[0] < 2 or values_b.shape[0] < 2:
                        continue

                    pval = ttest_ind(
                        values_a.values,
                        values_b.values,
                    ).pvalue

                    if pval < 0.05:
                        plot_ax.annotate(
                            '',
                            xy=(
                                index_a,
                                y_max + offset,
                            ),
                            xytext=(
                                index_b,
                                y_max + offset,
                            ),
                            xycoords='data',
                            textcoords='data',
                            arrowprops=dict(
                                arrowstyle='-',
                                ec='#000000',
                            ),
                        )
                        plot_ax.text(
                            x=np.mean([index_a, index_b]),
                            y=y_max + offset + y_max / 40,
                            s=pyp.utils.stars(pval, ns='-'),
                            horizontalalignment='center',
                            verticalalignment='center',
                        )
                        offset += y_max / 10

                plot_ax.set_ylim(
                    bottom=plot_ax.get_ylim()[0],
                    top=max([
                        plot_ax.get_ylim()[1],
                        y_max + offset + y_max / 10,
                    ]),
                )

        plot_ax.set_xlabel('')
        plot_ax.set_ylabel(
            '{} Signal'.format(
                'Relative' if cmp_groups else 'Cumulative',
            ),
        )
        plot_ax.get_legend().set_visible(False)

        if log_2:
            plot_ax.set_yticklabels(
                [
                    '{:.2f}'.format(i)
                    for i in np.power(2, plot_ax.get_yticks())
                ],
            )

        plot_ax.set_xticklabels(
            labels,
            rotation=45,
            horizontalalignment='right',
        )

        figures.append((plot_ax.get_figure(), plot_ax))

    return figures


def plot_all(
    data,
    cmp_groups=None,
):
    '''
    Runs :func:`.plot` and :func:`.plot_group` for all peptides in a data set.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.DataSet`
    cmp_groups : list of tuple, optional

    Returns
    -------
    figs : list of :class:`matplotlib.figure.Figure`
    '''
    figures = []

    figures += plot(
        data,
    )

    figures += plot_group(
        data,
        cmp_groups=cmp_groups,
    )

    return figures
