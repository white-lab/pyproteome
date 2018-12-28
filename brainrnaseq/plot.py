from collections import OrderedDict
import logging
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
import seaborn as sns

import brainrnaseq as brs

LOGGER = logging.getLogger('brainrnaseq.plot')


def plot_cell_enrichments(ds, enrichments, f=None, ax=None):
    LOGGER.info("Plotting cell type enrichments")

    if f is None:
        f = {'p': .05, 'asym_fold': 1.25}

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
    for cell in brs.CELL_TYPES:
        d = ds.filter(
            protein=[j for i in cell_prots.values() for j in i],
            fn=lambda x: len(x['Proteins']) < 2,
        )

        dc = d.filter(protein=cell_prots[cell])
        fore_hits = dc.filter(f).shape[0]
        fore_size = dc.shape[0]

        back_hits = d.filter(f).shape[0]
        back_size = d.shape[0]

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
                    ('p-value', -np.log10(val)),
                    ('hue', brs.CELL_COLORS[cell]),
                ])
            )
        )

    df = pd.DataFrame(vals)
    df = df.sort_values('p-value', ascending=False)

    ax = sns.barplot(
        data=df,
        y='cell',
        x='p-value',
        palette=df['hue'],
        ax=ax,
    )
    ax.set_title('Cell Type Enrichment')
    ax.set_ylabel('')
    ax.set_xlabel('p-value')
    ax.set_xticklabels(['{:.3}'.format(10 ** -i) for i in ax.get_xticks()])
    ax.axvline(-np.log10(.01), color='k', linestyle=':')

    return ax.get_figure(), ax
