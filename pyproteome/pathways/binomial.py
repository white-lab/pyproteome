
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


def _get_mods(ds, mods, accessions=False):
    for ind, row in ds.psms.iterrows():
        for prot in row['Proteins']:
            for mod in row['Modifications'].get_mods(mods):
                acc = prot.accession if accessions else prot.gene

                mod = '{}{}'.format(
                    mod.letter,
                    mod.abs_pos[
                        [
                            i.protein
                            for i in mod.sequence.protein_matches
                        ].index(prot)
                    ] + 1,
                )

                yield acc, mod

def _get_set(ds, mods=None, accessions=False):
    if mods:
        ret = _get_mods(ds, mods=mods, accessions=accessions)
    else:
        ret = ds.accessions if accessions else ds.genes

    return set(ret)


def plot_binomial_enrichment(
    ds, 
    gene_sets,
    filters=None,
    correlates=None, 
    corr_cutoff=.5, 
    fold_cutoff=1.5, 
    p_cutoff=1e-2,
    gene_rows=True,
    accessions=False,
    mods=None,
    minmax=4,
    **kwargs
):
    '''
    Plot binomial enrichment scores of arbitrary gene set lists.
    '''
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
            back = _get_set(
                fore, 
                accessions=accessions, 
                mods=mods,
            )

            if p_cutoff:
                fore = fore.filter(p=p_cutoff)

            if corr_cutoff:
                fore_up_ds = fore.filter(fn=lambda x: x['Fold Change'] > corr_cutoff)
                fore_down_ds = fore.filter(fn=lambda x: x['Fold Change'] < -corr_cutoff)
            else:
                fore_up_ds = fore.filter(asym_fold=fold_cutoff)
                fore_down_ds = fore.filter(asym_fold=1/fold_cutoff)

            fore_up = _get_set(
                fore_up_ds,
                accessions=accessions, 
                mods=mods,
            )
            fore_down = _get_set(
                fore_down_ds,
                accessions=accessions, 
                mods=mods,
            )

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
    
    cp = df.copy()
    cp['Group'] = cp.apply(lambda x: x['Group'] + '-' + x['Direction'] + '-' + x['Subset'], axis=1)
    del cp['Direction']
    del cp['Subset']
    
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

    cp.columns = cp.columns.droplevel()

    gene_colors = []
    n_sub = len(subsets)

    if correlates and len(correlates) > 1:
        colors = [
            'cyan',
            'magenta',
            'yellow',
            'blue',
        ]
        colors = sns.color_palette(
            'hls', len(correlates),
        ).as_hex()
        gene_colors += [[
            colors[i]
            for i in range(len(correlates))
            for j in range(n_sub)
        ]]
    
    gene_colors += [
        # (
        #     ['blue'] * len(subsets) + ['yellow'] * len(subsets)
        # ) * len(correlates),
        ['orange', 'green', 'red'][:n_sub] * (
            # 2 * len(correlates)
            len(correlates)
        ),
    ]

    figsize = (
        len(cp.columns) * .75,
        len(cp.index) / n_sub * 1.5,
    )

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

    cp = cp.loc[
        sorted(
            cp.index,
            key=lambda x: (
                list(correlates.keys()).index(x.split('-')[0]),
                list(filters.keys()).index(x.split('-')[-1]),
            )
        )
    ]
    
    if gene_rows:
        cp = cp.T
        cluster_kwargs['col_colors'] = gene_colors
    else:
        cluster_kwargs['row_colors'] = gene_colors

    display(cp)
    print(cluster_kwargs['col_colors'])

    g = sns.clustermap(
        cp,
        cbar_kws=dict(
            orientation='vertical',
            ticks=[-minmax, -minmax/2, 0, minmax/2, minmax],
        ),
        **cluster_kwargs
    )

    # g.cax.set_position([.75, .6, .02, .15])
    g.cax.set_title('LOE', loc='left', pad=10)

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
        (ax.axvline if gene_rows else ax.axhline)(v, ls='-', color='#4C4D4F', lw=1)
    
    return g, df
