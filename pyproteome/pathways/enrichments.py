# -*- coding: utf-8 -*-
'''
This module does most of the heavy lifting for the pathways module.

It includes functions for calculating enrichment scores and generating plots
for GSEA / PSEA.
'''
from __future__ import division

from collections import defaultdict
from functools import partial
import logging
import multiprocessing

import numpy as np
import pandas as pd
from sklearn.utils import shuffle

from . import plot


LOGGER = logging.getLogger('pyproteome.enrichments')
MIN_PERIODS = 5
'''
Minimum number of samples with peptide quantification and phenotypic
measurements needed to generate a correlation metric score.
'''
CORRELATION_METRICS = [
    'spearman',
    'pearson',
    'kendall',
    'fold',
    'log2',
    'zscore',
]
'''
Correlation metrics used for enrichment analysis. 'spearman', 'pearson', and
'kendall' are all calculated using `pandas.Series.corr()`.

'fold' takes ranking values direction from the 'Fold Change' column.

'log2' takes ranking values from a log2 'Fold Change' column.

'zscore' takes ranking values from a log2 z-scored 'Fold Change' column.
'''
DEFAULT_P = 0.75
DEFAULT_ESS_METHOD = 'integral'
DEFAULT_Q = 0.25
DEFAULT_RANK_CPUS = 6
'''
Default number of CPUs to use when scrambling rows of a data set.
'''
DEFAULT_CORR_CPUS = 4
'''
Default number of CPUs to use when scrambling columns of a data set.
'''
ESS_METHODS = {
    'max_abs': lambda x: max(x, key=abs),
    'max_min': lambda x: (max(x) - min(x)) * np.sign(max(x, key=abs)),
    'integral': lambda x: np.trapz(x),
}


class PrPDF(object):
    '''
    An exact probability distribution estimator.
    '''
    def __init__(self, data):
        self.data = np.sort(data)

    def pdf(self, x):
        '''
        Probability density function.

        Parameters
        ----------
        x : float

        Returns
        -------
        float
        '''
        return 1 / self.data.shape[0]

    def cdf(self, x):
        '''
        Cumulative density function.

        Parameters
        ----------
        x : float

        Returns
        -------
        float
        '''
        return np.searchsorted(self.data, x, side='left') / self.data.shape[0]

    def sf(self, x):
        '''
        Survival function.

        Parameters
        ----------
        x : float

        Returns
        -------
        float
        '''
        return 1 - (
            np.searchsorted(self.data, x, side='right') / self.data.shape[0]
        )


def _calc_p(ess, ess_pi):
    return [
        abs(ess) <= abs(i)
        for i in ess_pi
        if (i < 0) == (ess < 0)
    ]


def _shuffle(ser):
    ind = ser.index
    ser = shuffle(ser)
    ser.index = ind
    return ser


def simulate_es_s_pi(
    vals,
    psms,
    gene_sets,
    phenotype=None,
    metric='spearman',
    p_iter=1000,
    n_cpus=None,
    **kwargs
):
    '''
    Simulate ES(S, pi) by scrambling the phenotype / correlation values for a
    data set and recalculating gene set enrichment scores.

    Parameters
    ----------
    vals : :class:`pandas.DataFrame`
    psms : :class:`pandas.DataFrame`
    gene_sets : :class:`pandas.DataFrame`
    phenotype : :class:`pandas.Series`, optional
    metric : str, optional
    p_iter : int, optional
    n_cpus : int, optional

    Returns
    -------
    df : :class:`pandas.DataFrame`
    '''
    assert metric in CORRELATION_METRICS

    LOGGER.info(
        'Calculating ES(S, pi) for {} gene sets'.format(len(gene_sets))
    )

    if n_cpus is None:
        n_cpus = DEFAULT_RANK_CPUS

        if metric in ['spearman', 'pearson', 'kendall']:
            n_cpus = DEFAULT_CORR_CPUS

    if n_cpus > 1:
        pool = multiprocessing.Pool(
            processes=n_cpus,
        )
        gen = pool.imap_unordered(
            partial(
                _calc_essdist,
                psms=psms.copy(),
                gene_sets=gene_sets,
                metric=metric,
                **kwargs
            ),
            [phenotype for _ in range(p_iter)],
        )
    else:
        gen = (
            _calc_essdist(
                phen=phenotype,
                psms=psms,
                gene_sets=gene_sets,
                metric=metric,
                **kwargs
            )
            for _ in range(p_iter)
        )

    ess_dist = defaultdict(list)

    for ind, ess in enumerate(
        gen,
        start=1,
    ):
        for key, val in ess.items():
            ess_dist[key].append(val)

        if ind % (p_iter // min([p_iter, 10])) == 0:
            LOGGER.info(
                '-- Calculated {}/{} pvals'.format(ind, p_iter)
            )

    LOGGER.info('Calculating ES(S, pi) using {} cpus'.format(n_cpus))

    vals['ES(S, pi)'] = vals.index.map(
        lambda row: ess_dist[row],
    )

    return vals


def _calc_q(nes, nes_pdf, nes_pi_pdf):
    if nes > 0:
        return (
            nes_pi_pdf.pdf(nes) + nes_pi_pdf.sf(nes)
        ) / (
            nes_pdf.pdf(nes) + nes_pdf.sf(nes)
        )
    else:
        return (
            nes_pi_pdf.pdf(nes) + nes_pi_pdf.cdf(nes)
        ) / (
            nes_pdf.pdf(nes) + nes_pdf.cdf(nes)
        )


def estimate_pq(vals):
    '''
    Estimate p- and q-values for an enrichment analysis using the ES(S, pi)
    values generated by `simulate_es_s_pi`.

    Parameters
    ----------
    vals : :class:`pandas.DataFrame`
    '''
    assert 'ES(S)' in vals.columns
    assert 'ES(S, pi)' in vals.columns

    vals = vals.copy()

    mask = vals['ES(S)'] > 0

    pos_ess = vals['ES(S)'][mask]
    neg_ess = vals['ES(S)'][~mask]

    pos_pi = vals['ES(S, pi)'].apply(np.array).apply(lambda x: x[x > 0])
    neg_pi = vals['ES(S, pi)'].apply(np.array).apply(lambda x: x[x < 0])

    pos_mean = pos_pi.apply(np.mean).abs()
    neg_mean = neg_pi.apply(np.mean).abs()

    pos_nes = pos_ess / pos_mean
    neg_nes = neg_ess / neg_mean

    # assert (pos_nes.isnull() | neg_nes.isnull()).all()
    # assert ((~pos_nes.isnull()) | (~neg_nes.isnull())).all()

    pos_pi_nes = pos_pi / pos_mean
    neg_pi_nes = neg_pi / neg_mean

    vals['NES(S)'] = pos_nes.fillna(neg_nes)
    vals['pos NES(S, pi)'] = pos_pi_nes
    vals['neg NES(S, pi)'] = neg_pi_nes

    pos_mat = (
        np.concatenate(pos_pi_nes.values)
        if pos_pi_nes.shape[0] > 0 else
        np.array([])
    )
    neg_mat = (
        np.concatenate(neg_pi_nes.values)
        if neg_pi_nes.shape[0] > 0 else
        np.array([])
    )

    plot.plot_nes_dist(
        vals['NES(S)'].values,
        np.concatenate([pos_mat, neg_mat]),
    )

    pos_pdf = PrPDF(pos_nes.dropna().values)
    neg_pdf = PrPDF(neg_nes.dropna().values)

    pos_pi_pdf = PrPDF(pos_mat)
    neg_pi_pdf = PrPDF(neg_mat)

    LOGGER.info('Generated NES(S) distributions')

    if vals.shape[0] > 0:
        vals['p-value'] = vals.apply(
            lambda x:
            _frac_true(
                _calc_p(x['ES(S)'], x['ES(S, pi)'])
            ),
            axis=1,
        )

        vals['q-value'] = vals.apply(
            lambda x:
            _calc_q(
                x['NES(S)'],
                pos_pdf if x['NES(S)'] > 0 else neg_pdf,
                pos_pi_pdf if x['NES(S)'] > 0 else neg_pi_pdf,
            ),
            axis=1,
        )

        LOGGER.info('Calculated p, q values')

    return vals


def _frac_true(x):
    return sum(x) / max([len(x), 1])


def get_gene_changes(psms):
    '''
    Extract the gene IDs and correlation values for each gene / phosphosite in
    a data set. Merge together duplicate IDs by calculating their mean
    correlation.

    Parameters
    ----------
    psms : :class:`pandas.DataFrame`
    '''
    LOGGER.info('Getting gene correlations ({} IDs)'.format(psms.shape[0]))

    gene_changes = psms[['ID', 'Correlation']].copy()

    if gene_changes.shape[0] > 0:
        gene_changes = gene_changes.groupby(
            'ID',
            as_index=False,
        ).agg({
            'Correlation': np.mean,
        })

    gene_changes = gene_changes.sort_values(by='Correlation', ascending=False)
    gene_changes = gene_changes.set_index(keys='ID')

    return gene_changes


def _calc_essdist(
    phen=None, 
    psms=None, 
    gene_sets=None, 
    metric='spearman',
    **kwargs
):
    assert metric in CORRELATION_METRICS

    if phen is not None:
        phen = _shuffle(phen)

    if metric in ['fold', 'zscore']:
        psms = psms.copy()
        psms['Fold Change'] = _shuffle(psms['Fold Change'])

    vals = enrichment_scores(
        psms,
        gene_sets,
        phenotype=phen,
        metric=metric,
        recorrelate=True,
        pval=False,
        **kwargs
    )

    return vals['ES(S)']


def correlate_phenotype(psms, phenotype=None, metric='spearman'):
    '''
    Calculate the correlation values for each gene / phosphosite in a data set.

    Parameters
    ----------
    psms : :class:`pandas.DataFrame`
    phenotype : :class:`pandas.Series`, optional
    metric : str, optional
        The correlation function to use. See CORRELATION_METRICS for a full
        list of choices.
    
    Returns
    -------
    psms : :class:`pandas.DataFrame`
    '''
    assert metric in CORRELATION_METRICS

    psms = psms.copy()

    if metric in ['spearman', 'pearson', 'kendall']:
        LOGGER.info(
            'Calculating correlations using metric \'{}\' (samples: {})'
            .format(metric, list(phenotype.index))
        )
        psms['Correlation'] = psms.apply(
            lambda row:
            phenotype.corr(
                pd.to_numeric(row[phenotype.index]),
                method=metric,
                min_periods=MIN_PERIODS,
            ),
            axis=1,
        )
    else:
        LOGGER.info(
            'Calculating ranks'
        )
        new = psms['Fold Change']

        if (
            metric in ['log2'] or
            (metric in ['zscore'] and (new > 0).all())
        ):
            new = new.apply(np.log2)

        if metric in ['zscore']:
            new = (new - new.mean()) / new.std()

        psms['Correlation'] = new

    return psms


def calculate_es_s(gene_changes, gene_set, p=None, n_h=None, ess_method=None):
    '''
    Calculate the enrichment score for an individual gene set.

    Parameters
    ----------
    gene_changes : :class:`pandas.DataFrame`
    gene_set : set of str
    p : float, optional
    n_h : int, optional
    ess_method : str, optional
        One of {'integral', 'max_abs', 'max_min'}.

    Returns
    -------
    dict
    '''
    if p is None:
        p = DEFAULT_P

    if ess_method is None:
        ess_method = DEFAULT_ESS_METHOD

    assert ess_method in ESS_METHODS

    n = len(gene_changes)

    gene_set = set(
        gene
        for gene in gene_set
        if gene in gene_changes.index
    )
    hits = gene_changes.index.isin(gene_set)
    hit_list = gene_changes[hits].index.tolist()

    if n_h is None:
        n_h = len(gene_set)

    n_r = gene_changes[hits]['Correlation'].apply(
        lambda x: abs(x) ** p
    ).sum()

    scores = [0] + [
        ((abs(val) ** p) / n_r)
        if hit else
        (-1 / (n - n_h))
        for hit, val in zip(hits, gene_changes['Correlation'])
    ]
    cumsum = np.cumsum(scores)

    ess = ESS_METHODS.get(ess_method)(cumsum)

    return {
        'hits': hits,
        'cumscore': cumsum,
        'ess': ess,
        'hit_list': hit_list,
    }


def calculate_es_s_ud(gene_changes, up_set, down_set, **kwargs):
    '''
    Calculate the enrichment score for an individual gene set.

    Parameters
    ----------
    gene_changes : :class:`pandas.DataFrame`
    gene_set : set of str
    kwargs : dict, optional
        See extra arguments passed to calculate_es_s.

    Returns
    -------
    dict
    '''
    up_set = set(
        gene
        for gene in up_set
        if gene in gene_changes.index
    )
    down_set = set(
        gene
        for gene in down_set
        if gene in gene_changes.index
    )

    vals = []

    for gene_set in [up_set, down_set]:
        vals.append(
            calculate_es_s(
                gene_changes,
                gene_set,
                n_h=len(up_set) + len(down_set),
                **kwargs
            )
        )

    up_hits = vals[0]['hits']
    down_hits = vals[1]['hits']
    hit_list = vals[0]['hit_list'] + vals[1]['hit_list']

    upcumsum = vals[0]['cumscore']
    downcumsum = vals[1]['cumscore']

    ess = (
        vals[0]['ess'] if vals[0]['hits'].any() else 0
    ) - (
        vals[1]['ess'] if vals[1]['hits'].any() else 0
    )

    return {
        'hits': up_hits,
        'down_hits': down_hits,
        'cumscore': upcumsum,
        'down_cumscore': downcumsum,
        'ess': ess,
        'hit_list': hit_list,
    }


def _get_set_cols(cols):
    for set_cols in [
        ['up_set', 'down_set'],
        ['set'],
    ]:
        if any([i in cols for i in set_cols]):
            return set_cols


def enrichment_scores(
    psms,
    gene_sets,
    pval=True,
    recorrelate=False,
    metric=None,
    phenotype=None,
    **kwargs
):
    '''
    Calculate enrichment scores for each gene set.

    p-values and q-values are calculated by scrambling the phenotypes assigned
    to each sample or scrambling peptides' fold changes, depending on the
    correlation metric used.

    Parameters
    ----------
    psms : :class:`pandas.DataFrame`
    gene_sets : :class:`pandas.DataFrame`
    pval : bool, optional
    recorrelate : bool, optional
    metric : str, optional
    phenotype : :class:`pandas.Series`, optional
    kwargs : dict, optional
        See extra arguments passed to `calculate_es_s` and `simulate_es_s_pi`.

    Returns
    -------
    df : :class:`pandas.DataFrame`
    '''
    if metric is None:
        if phenotype is None:
            metric = 'fold'
        else:
            metric = 'spearman'

    assert metric in CORRELATION_METRICS

    if recorrelate:
        psms = correlate_phenotype(psms, phenotype=phenotype, metric=metric)

    gene_changes = get_gene_changes(psms)

    set_cols = _get_set_cols(gene_sets.columns)

    cols = [
        'name',
        'cumscore',
        'down_cumscore',
        'ES(S)',
        'hits',
        'down_hits',
        'hit_list',
        'n_hits',
    ] + set_cols
    vals = pd.DataFrame(
        columns=cols,
    )

    ess_args = {
        i: kwargs[i] 
        for i in ['p', 'ess_method'] 
        if i in kwargs
    }

    for set_id, row in gene_sets.iterrows():
        if set(set_cols) == set(['set']):
            es_vals = calculate_es_s(
                gene_changes,
                row['set'],
                **ess_args
            )
        else:
            es_vals = calculate_es_s_ud(
                gene_changes,
                row['up_set'],
                row['down_set'],
                **ess_args
            )

        vals = vals.append(
            pd.Series(
                [
                    row['name'],
                    es_vals.get('cumscore', []),
                    es_vals.get('down_cumscore', []),
                    es_vals.get('ess', 0),
                    es_vals.get('hits', []),
                    es_vals.get('down_hits', []),
                    es_vals.get('hit_list', []),
                    len(es_vals.get('hit_list', [])),
                ] + [row[i] for i in set_cols],
                name=set_id,
                index=cols,
            )
        )

    ess_pi_args = {
        i: kwargs[i] 
        for i in ['p', 'p_iter', 'n_cpus', 'ess_method'] 
        if i in kwargs
    }

    if pval:
        vals = simulate_es_s_pi(
            vals, psms, gene_sets,
            phenotype=phenotype,
            metric=metric,
            **ess_pi_args
        )

        vals = estimate_pq(vals)
        vals = vals.sort_values('NES(S)', ascending=False)
    else:
        vals = vals.sort_values('ES(S)', ascending=False)

    return vals


def filter_gene_sets(gene_sets, psms, min_hits=10):
    '''
    Filter gene sets to include only those with at least a given number of
    hits in a data set.

    Parameters
    ----------
    gene_sets : :class:`pandas.DataFrame`
    psms : :class:`pandas.DataFrame`
    min_hits : int, optional

    Returns
    -------
    df : :class:`pandas.DataFrame`
    '''
    LOGGER.info('Filtering gene sets')

    total_sets = gene_sets
    all_genes = set(psms['ID'])

    set_cols = _get_set_cols(gene_sets.columns)
    combined_set = gene_sets.apply(
        lambda x:
        set(
            i
            for col in set_cols
            for i in x[col]
        ),
        axis=1,
    )

    gene_sets = gene_sets[
        combined_set.apply(
            lambda x:
            len(x) < len(all_genes) and
            len(all_genes.intersection(x)) >= min_hits
        )
    ]

    LOGGER.info(
        'Filtered {} gene sets down to {} with ≥ {} genes present'
        .format(total_sets.shape[0], gene_sets.shape[0], min_hits)
    )

    return gene_sets


def filter_vals(
    vals,
    min_hits=0,
    min_abs_score=0,
    max_pval=1,
    max_qval=1,
):
    '''
    Filter gene set enrichment scores using give ES(S) / p-value / q-value
    cutoffs.

    Parameters
    ----------
    vals : :class:`pandas.DataFrame`
    min_hits : int, optional
    min_abs_score : float, optional
    max_pval : float, optional
    max_qval : float, optional

    Returns
    -------
    df : :class:`pandas.DataFrame`
    '''
    n_len = len(vals)
    if vals.shape[0] > 0:
        vals = vals[
            vals.apply(
                lambda x:
                x['n_hits'] >= min_hits and
                abs(x['ES(S)']) >= min_abs_score and (
                    (x['p-value'] <= max_pval)
                    if 'p-value' in x.index else
                    True
                ) and (
                    (x['q-value'] <= max_qval)
                    if 'q-value' in x.index else
                    True
                ),
                axis=1,
            )
        ]

    LOGGER.info(
        'Filtered {} gene sets down to {} after cutoffs'
        .format(n_len, len(vals))
    )

    return vals
