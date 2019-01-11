from collections import OrderedDict
import logging
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
import seaborn as sns

from pyproteome.motifs.plogo import format_title

import brainrnaseq as brs

LOGGER = logging.getLogger('brainrnaseq.plot')


def plot_cell_enrichments(
    ds,
    f=None,
    enrichments=None,
    ax=None,
    title=None,
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

    if ax is None:
        _, ax = plt.subplots(figsize=(4, 3))

    cell_prots = {
        cell: [key for key, val in enrichments.items() if val[0] == cell]
        for cell in brs.CELL_TYPES
    }

    display_name = {
        'Myelinating Oligodendrocytes': 'Oligodendrocytes',
    } if sum([
        'Oligo' in cell and len(cell_prots.get(cell, [])) > 0
        for cell in brs.CELL_TYPES
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

    for cell in brs.CELL_TYPES:
        for ind, fil in enumerate(f):
            dc = ds.filter(protein=set(cell_prots[cell]))

            fore_hits = dc.filter(fil).shape[0]
            fore_size = dc.shape[0]

            back_hits = ds.filter(fil).shape[0]
            back_size = ds.shape[0]

            if fore_size < 1 or back_size < 1:
                continue

            val = (
                1 -
                hypergeom.cdf(fore_hits, back_size, back_hits, fore_size) +
                hypergeom.pmf(fore_hits, back_size, back_hits, fore_size)
            )
            vals.append(
                pd.Series(
                    OrderedDict([
                        ('cell', display_name.get(cell, cell)),
                        ('fore hits', fore_hits),
                        ('fore size', fore_size),
                        ('back hits', back_hits),
                        ('back size', back_size),
                        ('p-value', val),
                        ('-log10 p-value', -np.log10(val)),
                        ('color', brs.CELL_COLORS[cell]),
                        ('hatch', hatches[ind % len(hatches)]),
                        ('hue', format_title(f=fil)),
                    ])
                )
            )

    df = pd.DataFrame(vals)

    ax = sns.barplot(
        data=df,
        y='cell',
        x='-log10 p-value',
        hue='hue',
        ax=ax,
    )

    ax.axvline(-np.log10(.01), color='k', linestyle=':')
    ax.legend(
        handles=[
            mpatches.Patch(
                facecolor='w',
                edgecolor='k',
                hatch=i,
                label=df['hue'].iloc[ind],
            )
            for ind, i in enumerate(hatches[:len(f)])
        ]
    )

    for hatch, color, p in zip(
        df['hatch'],
        df['color'],
        sorted(ax.patches, key=lambda x: x.xy[1]),
    ):
        p.set_hatch(hatch)
        p.set_facecolor(color)
        p.set_edgecolor('k')

    if title:
        ax.set_title(title)

    ax.set_ylabel('')
    ax.set_xlabel('p-value')
    ax.set_xticklabels(['{:.3}'.format(10 ** -i) for i in ax.get_xticks()])

    return ax.get_figure(), ax
