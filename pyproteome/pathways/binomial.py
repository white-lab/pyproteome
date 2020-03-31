
import numpy as np
from scipy import stats
import seaborn as sns
import pandas as pd

from pyproteome import data_sets


def log_odds(sf, cdf):
    if cdf == 0:
        score = -np.inf
    elif sf == 0:
        score = np.inf
    else:
        score = -np.log10(sf / cdf)
    
    return score


def binomial_scores(fore_up, fore_down, back, gene_sets):
    for key, vals in gene_sets.items():
        # Count of genes in each subset that exist in a gene set
        ct_back = [i for i in back if i in vals]
        ct_fore_up = [i for i in fore_up if i in vals]
        ct_fore_down = [i for i in fore_down if i in vals]

        # No hits
        if len(back) < 1 or len(ct_back) < 1:
            yield key, np.nan, np.nan
            continue

        # Calculate hypergeometric survival / cumulative distribution functions
        up_sf = stats.hypergeom.sf(
            len(ct_fore_up),
            len(back),
            len(ct_back),
            len(fore_up),
        )
        up_cdf = stats.hypergeom.cdf(
            len(ct_fore_up) - 1,
            len(back),
            len(ct_back),
            len(fore_up),
        )
        down_sf = stats.hypergeom.sf(
            len(ct_fore_down),
            len(back),
            len(ct_back),
            len(fore_down),
        )
        down_cdf = stats.hypergeom.cdf(
            len(ct_fore_down) - 1,
            len(back),
            len(ct_back),
            len(fore_down),
        )
        # Scores are log odds-ratio values: -log10(P(foreground) / P(background))
        up_score = log_odds(up_sf, up_cdf)
        down_score = log_odds(down_sf, down_cdf)
        if down_cdf == 0:
            down_score = -np.inf
        elif down_sf == 0:
            down_score = np.inf
        else:
            down_score = -np.log10(down_sf / down_cdf)

        yield key, up_score, down_score


def plot_binomial_enrichment(
    ds, 
    gene_sets,
    filters=None,
    correlates=None, 
    corr_cutoff=.5, 
    fold_cutoff=1.5, 
    p_cutoff=1e-2,
    gene_rows=True,
    minmax=4,
    **kwargs
):
    """
    Plot binomial enrichment scores of arbitrary gene set lists.
    """
    if correlates is None:
        corr_cutoff = None
        correlates = {'dummy': None}

    df = []

    subsets = {
        name:
        data_sets.merge_proteins(ds.filter(f)) if merge else ds.filter(f)
        for name, (f, merge) in filters.items()
    }

    for axis_name, target in correlates.items():
        for sub_name, sub_cp in subsets.items():
            fore = sub_cp.copy()

            if corr_cutoff:
                fore = data_sets.update_correlation(fore, target, metric='spearman')

            fore.psms = fore.psms.dropna(subset=('Fold Change', 'p-value'))
            back = set(fore.genes)

            if p_cutoff:
                fore = fore.filter(p=p_cutoff)

            if corr_cutoff:
                fore_up = set(fore.filter(fn=lambda x: x['Fold Change'] > corr_cutoff).genes)
                fore_down = set(fore.filter(fn=lambda x: x['Fold Change'] < -corr_cutoff).genes)
            else:
                fore_up = set(fore.filter(asym_fold=fold_cutoff).genes)
                fore_down = set(fore.filter(asym_fold=1/fold_cutoff).genes)

            for set_name, up_score, down_score in binomial_scores(
                fore_up, fore_down, back, gene_sets,
            ):
                df.append((
                    axis_name,
                    set_name,
                    sub_name,
                    'correlated' if corr_cutoff else 'upregulated',
                    up_score,
                ))
                df.append((
                    axis_name,
                    set_name,
                    sub_name,
                    'anticorrelated' if corr_cutoff else 'downregulated',
                    down_score,
                ))
            
    df = pd.DataFrame(
        df,
        columns=['Group', 'Gene Set', 'Subset', 'Direction', 'Score'],
    )
    df = df[df['Direction'].isin(['correlated', 'upregulated'])]

    # f, ax = plt.subplots(
    #     dpi=200,
    #     figsize=(3.5, 12),
    # )
    
    cp = df.copy()
    cp['Group'] = cp.apply(lambda x: x['Group'] + '-' + x['Direction'] + '-' + x['Subset'], axis=1)
    del cp['Direction']
    del cp['Subset']

    display(cp[cp['Score'] == np.inf])
    
    cp = cp.pivot(
        index='Group', 
        columns='Gene Set',
    ).fillna(
        0,
    ).replace(
        [-np.inf], 0,
    ).replace(
        [np.inf], minmax,
        # [np.inf], 0,
    ).iloc[::-1]

    gene_colors = []

    if correlates and len(correlates) > 1:
        gene_colors += [
            ['magenta'] * (len(correlates) * len(subsets) // 2) + 
            ['cyan'] * (len(correlates) * len(subsets) // 2),
        ]
    
    gene_colors += [
        # (
        #     ['blue'] * len(subsets) + ['yellow'] * len(subsets)
        # ) * len(correlates),
        ['orange', 'green', 'red'] * (
            # 2 * len(correlates)
            len(correlates)
        ),
    ]
    print(list(correlates.keys()))
    print(list(subsets.keys()))
    print(gene_colors)

    figsize = (
        len(cp.index) * 2,
        len(cp.columns) / 10 / 2, 
    )
    # figsize = (
    #     len(cp.index) * 1.2,
    #     len(cp.columns),
    # )

    if gene_rows:
        figsize = tuple(reversed(figsize))

    cluster_kwargs = {
        'cmap': 'coolwarm',
        'vmin': -minmax,
        'vmax': minmax,
        'row_cluster': False,
        'col_cluster': False,
        'figsize': figsize,
    }
    cluster_kwargs.update(kwargs)
    
    if gene_rows:
        cp = cp.T
        cluster_kwargs['col_colors'] = gene_colors
    else:
        cluster_kwargs['row_colors'] = gene_colors

    g = sns.clustermap(
        cp,
        cbar_kws=dict(
            label='log odds enrichment',
    #         use_gridspec=False,
    #         location="top",
            orientation="vertical",
        ),
        **cluster_kwargs
    )

    g.cax.set_position([.92, .27, .03, .3])

    g.ax_row_dendrogram.set_visible(False)
    g.ax_col_dendrogram.set_visible(False)

    ax = g.ax_heatmap

    ax.tick_params(
        top=False, 
        bottom=False, 
        left=False, 
        right=False, 
        labelleft=True, 
        labelbottom=True,
        labelright=False,
        labeltop=False,
    )

    ax.set_xlabel('')
    ax.set_ylabel('')

    (ax.set_xticklabels if gene_rows else ax.set_yticklabels)([])
    (ax.set_yticklabels if gene_rows else ax.set_xticklabels)(
        [
            i.get_text().split('-', 1)[-1]#.split(' ')[-1]
            for i in (ax.get_yticklabels() if gene_rows else ax.get_xticklabels())
        ],
        rotation=0 if gene_rows else 45,
        ha='right',
    )
    # ax.get_colorbar()
    for v in range(
        len(subsets), 
        len(cp.columns if gene_rows else cp.index), 
        len(subsets)
    ):
        (ax.axvline if gene_rows else ax.axhline)(v, ls='-', color='k', lw=1)
    
    return g, df
