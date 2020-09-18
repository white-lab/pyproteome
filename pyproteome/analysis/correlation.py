# -*- coding: UTF-8 -*-
'''
This module provides functionality for data set analysis.

Functions include volcano plots, sorted tables, and plotting sequence levels.
'''

from __future__ import absolute_import, division

import logging
import warnings

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import pyproteome as pyp


LOGGER = logging.getLogger('pyproteome.correlation')


def correlate_data_sets(
    data1, data2,
    adjust=True,
    label_cutoff=1.5,
    show_labels=False,
    show_title=True,
    ax=None,
):
    '''
    Plot the correlation between peptides levels in two different data sets.

    Parameters
    ----------
    data1 : :class:`pyproteome.data_sets.DataSet`
    data2 : :class:`pyproteome.data_sets.DataSet`
    filename : str, optional
    '''
    LOGGER.info('Plotting dataset correlation')

    merged = pd.merge(
        data1.psms, data2.psms,
        on='Sequence',
    ).dropna(subset=('Fold Change_x', 'Fold Change_y'))

    if ax is None:
        f, ax = plt.subplots(
            figsize=(4, 3),
        )
    else:
        f = ax.get_figure()

    merged['Fold Change_x'] = merged['Fold Change_x'].apply(np.log2)
    merged['Fold Change_y'] = merged['Fold Change_y'].apply(np.log2)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=FutureWarning)
        sns.regplot(
            x='Fold Change_x',
            y='Fold Change_y',
            data=merged,
            ax=ax,
            scatter_kws={
                's': 2,
            },
        )

    label_cutoff = np.log2(label_cutoff)

    if show_labels:
        texts = []

        for _, row in merged.iterrows():
            x = row['Fold Change_x']
            y = row['Fold Change_y']
            ratio = x - y

            if ratio < label_cutoff and ratio > - label_cutoff:
                continue

            txt = ' / '.join(row['Proteins_x'].genes)
            txt = txt[:20] + ('...' if len(txt) > 20 else '')

            text = ax.text(
                x=x,
                y=y,
                s=txt,
                bbox=dict(
                    color='lightgreen' if ratio < 0 else 'pink',
                    alpha=0.8,
                ),
            )

            texts.append(text)

        if adjust and texts:
            pyp.utils.adjust_text(
                texts=texts,
                ax=ax,
                lim=400,
                force_text=0.1,
                force_points=0.1,
                arrowprops=dict(arrowstyle='->', relpos=(0, 0), lw=1),
                only_move={
                    'points': 'y',
                    'text': 'xy',
                }
            )

    min_x = merged['Fold Change_x'].min()
    max_x = merged['Fold Change_x'].max()
    min_y = merged['Fold Change_y'].min()
    max_y = merged['Fold Change_y'].max()

    min_xy = min([min_x, min_y])
    max_xy = max([max_x, max_y])

    ax.plot([min_xy, max_xy], [min_xy, max_xy], '--', color='k')
    ax.set_xlim(left=min_xy - .5, right=max_xy + .5)
    ax.set_ylim(bottom=min_xy - .5, top=max_xy + .5)

    name1 = data1.name
    name2 = data2.name

    ax.set_xlabel('Fold Change -- {}'.format(name1))
    ax.set_ylabel('Fold Change -- {}'.format(name2))

    ax.set_xticklabels(
        ['{:.2f}'.format(i) for i in np.power(2, ax.get_xticks())],
    )
    ax.set_yticklabels(
        ['{:.2f}'.format(i) for i in np.power(2, ax.get_yticks())],
    )

    if show_title:
        pear_corr = merged['Fold Change_x'].corr(
            merged['Fold Change_y'],
            method='pearson',
        )
        # spear_corr = merged['Fold Change_x'].corr(
        #     merged['Fold Change_y'],
        #     method='spearman',
        # )

        ax.set_title(
            (
                r'Pearson\'s: $\rho$={:.2f}'
                # r'Spearman\'s: $\rho$={:.2f}'
            ).format(pear_corr)
        )

    return f


def _scatter_plots(
    cp, signal,
    data_chans, signal_chans, signal_groups,
    scatter_cols=4,
    xlabel='',
    scatter_colors=None,
    scatter_symbols=None,
):
    scatter_colors = scatter_colors or {}
    scatter_symbols = scatter_symbols or {}

    if cp.shape[0] < 1:
        return None

    df = pd.DataFrame(
        [
            (
                row,
                str(row['Sequence']),
                signal[sig_chan],
                row[data_chan],
                sig_group,
                scatter_colors.get(sig_group, 'black'),
            )
            for _, row in cp.psms.sort_values('Correlation').iterrows()
            for data_chan, sig_chan, sig_group in zip(
                data_chans, signal_chans, signal_groups,
            )
        ],
        columns=('row', 'sequence', 'x', 'y', 'group', 'color'),
    )
    g = sns.lmplot(
        x='x',
        y='y',
        hue='group',
        col='sequence',
        col_wrap=scatter_cols,
        data=df,
        sharey=False,
        fit_reg=False,
    )
    for ax in g.axes:
        seq = ax.get_title().split('=', 1)[1].strip()
        df_cp = df[df['sequence'] == seq]
        row = df_cp['row'].iloc[0]

        row_title = ' / '.join(str(i.gene) for i in row['Proteins'])
        row_seq = str(row['Sequence'])
        row_mods = str(row['Modifications'].get_mods([(None, 'Phospho')]))

        sns.regplot(
            x='x',
            y='y',
            data=df_cp,
            scatter=False,
            ax=ax,
        )

        ax.set_title(
            row_title[:20] + ('...' if len(row_title) > 20 else ''),
            # fontsize=28,
            # fontweight='bold',
        )
        # ax.set_xlabel(
        #     '{}\n$\\rho = {:.2f}$'.format(
        #         xlabel,
        #         row['Correlation'],
        #     ),
        #     fontsize=22,
        # )
        ax.set_ylabel(
            row_seq[:20] + (
                '...' if len(row_seq) > 20 else ''
            ) + (
                '\n({})'.format(
                    row_mods[:20] + ('...' if len(row_mods) > 20 else '')
                ) if row_mods else ''
            ),
            # fontsize=20,
        )

        # for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
        #     tick.label.set_fontsize(20)

        ax.legend(loc='best')

    return g


def _remove_lesser_dups(labels, compress_sym=False):
    new_labels = []
    labels = list(labels)

    for index, (x, y, label) in enumerate(labels):
        for o_index, (o_x, o_y, o_label) in enumerate(labels):
            if index == o_index:
                continue
            if label != o_label:
                continue
            if not compress_sym and (y < 0) != (y < 0):
                continue

            mul = 1 if y >= 0 else -1
            o_mul = 1 if o_y >= 0 else -1

            if (
                mul * y < o_mul * o_y
            ) or (
                (
                    mul * y <= o_mul * o_y
                ) and
                index < o_index
            ):
                break
        else:
            new_labels.append((x, y, label))

    return new_labels


def correlate_signal(
    data,
    signal,
    corr_cutoff=0.8,
    scatter_cols=4,
    options=None,
    title=None,
    show_duplicates=False,
    scatter_colors=None,
    scatter_symbols=None,
    show_scatter=True,
    ax=None,
    xlabel='',
):
    '''
    Calculate the correlation between levels of each peptide in a data set and
    a given signal variable.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.DataSet`
    signal : :class:`pandas.Series`

    Returns
    -------
    f_corr : :class:`matplotlib.figure.Figure`
    f_scatter : :class:`matplotlib.figure.Figure`
    '''
    options = options or {}

    highlight = options.get('highlight', {})
    hide = options.get('hide', {})
    edgecolors = options.get('edgecolors', {})
    rename = options.get('rename', {})
    scatter_colors = options.get('scatter_colors', {})
    scatter_symbols = options.get('scatter_symbols', {})

    if title is None:
        title = data.name

    cp = data.copy()
    cp = cp.rename_channels()

    cols = getattr(signal, 'columns', getattr(signal, 'index', []))

    signal_groups = [
        label
        for label, group in cp.groups.items()
        for chan in group
        if chan in cp.channels.keys() and
        chan in cols
    ]

    signal_chans = [
        chan
        for group in cp.groups.values()
        for chan in group
        if chan in cp.channels.keys() and
        chan in cols
    ]
    data_chans = [
        cp.channels[chan]
        for chan in signal_chans
    ]

    cp.psms['Correlation'] = cp.psms.apply(
        lambda row:
        signal[signal_chans].corr(
            pd.to_numeric(row[data_chans]),
            method='spearman',
            min_periods=5,
        ),
        axis=1,
    ).replace([np.inf, -np.inf], np.nan)
    cp.psms = cp.psms.dropna(
        subset=['Correlation'],
        how='any',
    )

    if ax is None:
        f_corr, ax = plt.subplots(figsize=(12, 10))

    x, y, colors = [], [], []
    sig_x, sig_y, sig_labels = [], [], []

    for index, (_, row) in enumerate(cp.psms.iterrows()):
        x.append(index)
        y.append(row['Correlation'])

        sig = abs(row['Correlation']) >= corr_cutoff

        if sig:
            sig_x.append(index)
            sig_y.append(row['Correlation'])
            sig_labels.append(' / '.join(row['Proteins'].genes))

        colors.append(
            'blue'
            if sig else
            '{:.2f}'.format(
                max([len(row[data_chans].dropna()) / len(data_chans) - .25, 0])
            )
        )

    ax.scatter(x, y, c=colors)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    labels = zip(sig_x, sig_y, sig_labels)

    if not show_duplicates:
        labels = _remove_lesser_dups(labels)

    labels = sorted(labels, key=lambda x: x[1], reverse=True)[:100]
    LOGGER.info('Showing names for {} genes'.format(len(labels)))

    texts = []

    for xs, ys, txt in labels:
        if any(i.strip() in hide for i in txt.split('/')):
            continue

        txt = rename.get(txt, txt)

        texts.append(
            ax.text(
                x=xs,
                y=ys,
                s=txt[:20] + ('...' if len(txt) > 20 else ''),
                fontsize=20 if txt in highlight else 16,
                bbox=dict(
                    facecolor='lightgreen' if xs > 0 else 'pink',
                    alpha=1,
                    linewidth=0.5 if txt not in edgecolors else 3,
                    edgecolor=edgecolors.get(txt, 'black'),
                    boxstyle='round',
                ),
            )
        )

    if texts:
        pyp.utils.adjust_text(
            texts=texts,
            ax=ax,
            lim=400,
            force_text=0.3,
            force_points=0.01,
            arrowprops=dict(arrowstyle='->', relpos=(0, 0), lw=1),
            only_move={
                'points': 'y',
                'text': 'xy',
            }
        )

    ax.set_xlabel(
        'Index',
    )
    # ax.set_yticklabels(
    #     '{:.3}'.format(i)
    #     for i in np.power(1/10, ax.get_yticks())
    # )
    ax.set_ylabel(
        'Correlation',
    )

    if title:
        ax.set_title(title)

    cp.psms = cp.psms[cp.psms['Correlation'].apply(abs) >= corr_cutoff]

    f_scatter = None

    if show_scatter:
        f_scatter = _scatter_plots(
            cp, signal, data_chans, signal_chans, signal_groups,
            scatter_cols=scatter_cols,
            scatter_colors=scatter_colors,
            scatter_symbols=scatter_symbols,
        )

    return ax.get_figure(), f_scatter
