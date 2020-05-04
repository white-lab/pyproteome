from collections import OrderedDict
import logging
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
import seaborn as sns

from pyproteome.motifs.plogo import format_title

import pyproteome as pyp
import brainrnaseq as brs

LOGGER = logging.getLogger('brainrnaseq.plot')


def plot_cell_enrichments(
    ds,
    f=None,
    enrichments=None,
    title=None,
    show_cbar=True,
    **kwargs
):
    LOGGER.info("Plotting cell type enrichments")

    if enrichments is None:
        enrichments = brs.enrichments.get_enrichments(
            list(ds.species)[0],
        )

    if f is None:
        f = {'p': .05, 'asym_fold': 1.25}

    if isinstance(f, dict):
        f = [f]

    cell_types = [
        cell
        for cell in brs.CELL_TYPES
        if cell not in set(['Endothelia', 'OPC', 'New Oligodendrocytes'])
    ]

    cell_prots = {
        cell: [key for key, val in enrichments.items() if val == cell]
        for cell in cell_types
    }

    display_name = {
        'New Oligodendrocytes': 'Oligodendrocytes',
        'Myelinating Oligodendrocytes': 'Oligodendrocytes',
        'OPC': 'Oligodendrocytes',
    } if sum([
        'Oligo' in cell and len(cell_prots.get(cell, [])) > 0
        for cell in cell_types
    ]) else {}

    vals = []

    ds = ds.filter(
        protein=set(j for i in cell_prots.values() for j in i),
        fn=lambda x: len(x['Proteins']) < 2,
    )
    hatches = [
        "",
        "//",
        "o",
        "x",
        ".",
    ]

    groups = pyp.utils.flatten_list(
        ds.cmp_groups or [i for i in [ds.group_a, ds.group_b] if i]
    )

    for ind, fil in enumerate(f):
        for cell_ind, cell in enumerate(cell_types):
            LOGGER.info(
                'Calculating enrichment for {} - {} ({}/{})'.format(
                    cell,
                    format_title(f=fil),
                    1 + ind * len(cell_types) + cell_ind,
                    len(f) * len(cell_types),
                )
            )

            for group in pyp.utils.flatten_list(groups):
                for sample in ds.groups[group]:
                    if sample not in ds.channels:
                        continue

                    dc = ds.filter(protein=set(cell_prots[cell]))

                    chan = dc.channels[sample]
                    dc.psms['Fold Change'] = dc[chan]
                    dc.psms = dc.psms.dropna(subset=['Fold Change'], axis=0)
                    fore_hits = dc.filter(fil).shape[0]
                    fore_size = dc.shape[0]

                    back_hits = ds.filter(fil).shape[0]
                    back_size = ds.shape[0]

                    if fore_size < 1 or back_size < 1:
                        continue

                    # if fore_size < 10:
                    #     continue

                    val = (
                        1 -
                        hypergeom.cdf(
                            fore_hits, back_size, back_hits, fore_size
                        ) +
                        hypergeom.pmf(
                            fore_hits, back_size, back_hits, fore_size
                        )
                    )
                    log_val = (
                        -np.log10(val) if val > 0 else np.nan
                    )
                    vals.append(
                        pd.Series(
                            OrderedDict([
                                ('cell', display_name.get(cell, cell)),
                                ('sample', sample),
                                ('chan', chan),
                                ('fore hits', fore_hits),
                                ('fore size', fore_size),
                                ('back hits', back_hits),
                                ('back size', back_size),
                                ('p-value', val),
                                ('-log10 p-value', log_val),
                                ('color', brs.CELL_COLORS[cell]),
                                ('hatch', hatches[ind % len(hatches)]),
                                ('hue', format_title(f=fil)),
                            ])
                        )
                    )

    df = pd.DataFrame(vals)
    df['cell'] = df.apply(
        lambda x: '{} - {}'.format(x['cell'], x['hue']),
        axis=1,
    )
    df = df.set_index(['cell', 'sample'])
    ind = df.index
    df = df.unstack()['-log10 p-value']
    # df = df.unstack()['fore hits']
    df = df[ind.get_level_values(1).unique()]
    # df = df.replace(np.inf, np.nan).fillna(0)
    print(df.max().max())
    df = df.fillna(0)
    df = df.replace(np.inf, np.nan).fillna(df.max().max())
    ax = sns.heatmap(
        df,
        cmap='viridis',
        cbar=False,
        **kwargs
    )

    if show_cbar:
        cbar = ax.figure.colorbar(
            ax.collections[0],
            cax=ax.inset_axes([0.05, -.6, .9, .05]),
            orientation="horizontal",
        )
        from matplotlib import ticker
        cbar.locator = ticker.MaxNLocator(nbins=4)
        cbar.update_ticks()
        cbar.set_label('-log$_{10}$ p-value')
        cbar.ax.xaxis.set_label_position('top')

        # cbar.ax.set_xticklabels(
        #     ['{:.2e}'.format(10 ** -i) for i in cbar.ax.get_xticks()]
        # )
        cbar.outline.set_visible(False)

    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    ax.set_xlabel('')
    ax.set_ylabel('')

    for ind in range(len(cell_types) - 1):
        ax.axhline((ind + 1) * len(f), color='grey', linestyle='-')

    xticks = []
    x = 0
    for group in groups:
        chan_len = len([
            sample
            for sample in ds.groups[group]
            if sample in ds.channels
        ])

        if chan_len < 1:
            continue

        xticks.append(x + (chan_len / 2))
        x += chan_len

        if x < ax.get_xlim()[1]:
            ax.axvline(x, color='grey', linestyle=":")

    ax.set_xticks(xticks)
    ax.set_xticklabels(
        [
            i
            for fil in f
            for i in groups
            if any([sample in ds.channels for sample in ds.groups[i]])
        ],
        rotation=45,
        ha='right',
    )

    # ax = sns.barplot(
    #     data=df,
    #     y='cell',
    #     x='-log10 p-value',
    #     hue='hue',
    #     ax=ax,
    # )
    #
    # ax.axvline(-np.log10(.01), color='k', linestyle=':')
    # ax.legend(
    #     handles=[
    #         mpatches.Patch(
    #             facecolor='w',
    #             edgecolor='k',
    #             hatch=i,
    #             label=df['hue'].iloc[ind],
    #         )
    #         for ind, i in enumerate(hatches[:len(f)])
    #     ]
    # )
    #
    # for hatch, color, p in zip(
    #     df['hatch'],
    #     df['color'],
    #     sorted(ax.patches, key=lambda x: x.xy[1]),
    # ):
    #     p.set_hatch(hatch)
    #     p.set_facecolor(color)
    #     p.set_edgecolor('k')
    #
    # if title:
    #     ax.set_title(title)
    #
    # ax.set_ylabel('')
    # ax.set_xlabel('p-value')
    # ax.set_xticklabels(['{:.3}'.format(10 ** -i) for i in ax.get_xticks()])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

    return ax.get_figure(), ax
