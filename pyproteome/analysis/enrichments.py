
from __future__ import division

from collections import defaultdict
from functools import partial
import logging
import multiprocessing

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.utils import shuffle


LOGGER = logging.getLogger("pyproteome.enrichments")
MIN_PERIODS = 5


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
    ds, gene_sets, phenotype,
    p=1,
    p_iter=1000,
):
    LOGGER.info(
        "Calculating ES(S, pi) for {} gene sets".format(len(gene_sets))
    )
    pool = multiprocessing.Pool(
        processes=6,
    )

    ess_dist = defaultdict(list)

    for ind, ess in enumerate(
        pool.imap_unordered(
            partial(_calc_essdist, ds=ds, gene_sets=gene_sets, p=p),
            [_shuffle(phenotype) for _ in range(p_iter)],
        ),
        start=1,
    ):
        for key, val in ess.items():
            ess_dist[key].append(val)

        if ind % (p_iter // 10) == 0:
            LOGGER.info(
                "Calculated {}/{} pvals".format(ind, p_iter)
            )

    LOGGER.info("Calculating ES(S, pi)")

    vals["ES(S, pi)"] = vals.index.map(
        lambda row: ess_dist[row],
    )

    return vals


def _calc_q(nes, nes_pdf, nes_pi_pdf):
    return (
        nes_pi_pdf.pdf(nes) + nes_pi_pdf.sf(nes)
    ) / (
        nes_pdf.pdf(nes) + nes_pdf.sf(nes)
    )


class PrPDF(object):
    def __init__(self, data):
        self.data = np.sort(data)

    def pdf(self, x):
        return 1 / self.data.shape[0]

    def cdf(self, x):
        return np.searchsorted(self.data, x, side="left") / self.data.shape[0]

    def sf(self, x):
        return 1 - self.cdf(x) - self.pdf(x)


def estimate_pq(vals):
    vals = vals.copy()

    pos_ess = vals["ES(S)"][vals["ES(S)"] > 0]
    neg_ess = vals["ES(S)"][vals["ES(S)"] < 0]

    pos_pi = vals["ES(S, pi)"].apply(np.array).apply(lambda x: x[x > 0])
    neg_pi = vals["ES(S, pi)"].apply(np.array).apply(lambda x: x[x < 0])

    pos_mean = pos_pi.apply(np.mean)
    neg_mean = neg_pi.apply(np.mean)

    pos_nes = (pos_ess / pos_mean)
    neg_nes = (neg_ess / neg_mean)

    pos_pi_nes = pos_pi / pos_mean
    neg_pi_nes = neg_pi / neg_mean

    vals["NES(S)"] = pos_nes.fillna(neg_nes)
    vals["pos NES(S, pi)"] = pos_pi_nes
    vals["neg NES(S, pi)"] = neg_pi_nes

    for index, row in pos_nes.items():
        assert not np.isnan(row) or not np.isnan(neg_nes[index])
    for index, row in neg_nes.items():
        assert not np.isnan(row) or not np.isnan(pos_nes[index])

    LOGGER.info("Calculated pos, neg distributions")

    pos_pdf = PrPDF(pos_nes.dropna().as_matrix())
    neg_pdf = PrPDF(neg_nes.dropna().as_matrix())

    pos_pi_pdf = PrPDF(np.concatenate(pos_pi_nes.as_matrix()))
    neg_pi_pdf = PrPDF(np.concatenate(neg_pi_nes.as_matrix()))

    LOGGER.info("Normalized NES distributions")

    vals["p-value"] = vals.apply(
        lambda x:
        _frac_true(
            _calc_p(x["ES(S)"], x["ES(S, pi)"])
        ),
        axis=1,
    )

    vals["q-value"] = vals.apply(
        lambda x:
        _calc_q(
            x["NES(S)"],
            pos_pdf if x["ES(S)"] > 0 else neg_pdf,
            pos_pi_pdf if x["ES(S)"] > 0 else neg_pi_pdf,
        ),
        axis=1,
    )

    LOGGER.info("Calculated p, q values")
    return vals


def _frac_true(x):
    return sum(x) / max([len(x), 1])


def absmax(i):
    return max([abs(min(i)), abs(max(i))])


def _get_changes(ds):
    gene_changes = ds.psms[["Entrez", "Correlation"]]
    gene_changes.is_copy = None

    gene_changes["Abs Corr"] = gene_changes["Correlation"].apply(abs)
    gene_changes = gene_changes.sort_values(by="Abs Corr", ascending=True)
    gene_changes = gene_changes.drop_duplicates(subset="Entrez", keep="last")

    gene_changes = gene_changes.sort_values(by="Correlation", ascending=True)
    gene_changes = gene_changes.set_index(keys="Entrez")

    return gene_changes


def _calc_essdist(phen, ds=None, gene_sets=None, p=1):
    vals = enrichment_scores(
        ds,
        gene_sets,
        phen,
        p=p,
        recorrelate=True,
        pval=False,
    )

    return vals["ES(S)"]


def correlate_phenotype(ds, phenotype):
    ds = ds.copy()
    ds.psms["Correlation"] = ds.psms.apply(
        lambda row:
        phenotype.corr(
            pd.to_numeric(row[phenotype.index]),
            method="pearson",
            min_periods=MIN_PERIODS,
        ),
        axis=1,
    )
    return ds


def calculate_es_s(gene_changes, gene_set, p=1):
    n = len(gene_changes)
    gene_set = [
        gene
        for gene in gene_set
        if gene in gene_changes.index
    ]
    hits = gene_changes.index.isin(gene_set)

    n_h = len(gene_set)
    n_r = gene_changes[hits]["Correlation"].apply(
        lambda x: abs(x) ** p
    ).sum()

    scores = [0] + [
        ((abs(val) ** p) / n_r)
        if hit else
        (-1 / (n - n_h))
        for hit, val in zip(hits, gene_changes["Correlation"])
    ]
    return hits, np.cumsum(scores)


def enrichment_scores(
    ds, gene_sets, phenotype,
    p=1, pval=True, recorrelate=False, p_iter=1000,
):
    """
    Here is the function for calculating gene set enrichment scores.
    It calculates the association of each gene with a given phenotype
    and then generates the ES(S) scores for a given gene set.
    """
    if recorrelate:
        ds = correlate_phenotype(ds, phenotype)

    gene_changes = _get_changes(ds)

    cols = ["cumscore", "ES(S)", "hits", "n_hits"]
    vals = pd.DataFrame(
        columns=cols,
    )

    for set_id, row in gene_sets.iterrows():
        hits, cumscore = calculate_es_s(
            gene_changes,
            row["set"],
            p=p,
        )
        vals = vals.append(
            pd.Series(
                [
                    cumscore,
                    max(cumscore, key=abs),
                    hits,
                    sum(hits),
                ],
                name=set_id,
                index=cols,
            )
        )

    if pval:
        vals = simulate_es_s_pi(
            vals, ds, gene_sets, phenotype,
            p=p,
            p_iter=p_iter,
        )
        vals = estimate_pq(vals)

    return vals


# Here's a helper function to generate the plots seen below
def plot_enrichment(
    ds, gene_sets, phenotype,
    p=1,
    cols=5,
    min_hits=10,
    min_abs_score=0.4,
    max_pval=1,
    max_qval=1,
    pval=False,
    p_iter=1000,
):
    # Plot the ranked list of correlations
    f, ax = plt.subplots()
    gene_changes = _get_changes(ds)

    ax.plot(sorted(
        gene_changes["Correlation"],
        reverse=True,
    ))
    ax.axhline(0, color="k")
    ax.set_xlabel("Gene List Rank", fontsize=20)
    ax.set_ylabel("Correlation", fontsize=20)

    gene_sets = gene_sets[
        (gene_sets["set"].apply(len) < len(set(ds["Entrez"]))) &
        gene_sets["set"].apply(
            lambda x: len(set(ds["Entrez"]).intersection(x)) > min_hits
        )
    ]

    vals = enrichment_scores(
        ds,
        gene_sets,
        phenotype,
        p=p,
        pval=pval,
        p_iter=p_iter,
    )

    filtered_vals = vals[
        vals.apply(
            lambda x:
            abs(x["ES(S)"]) >= min_abs_score and (
                (x["p-value"] < max_pval)
                if "p-value" in x.index else
                True
            ) and (
                (x["q-value"] < max_qval)
                if "q-value" in x.index else
                True
            ),
            axis=1,
        )
    ]
    filtered_vals = filtered_vals.sort_values("ES(S)", ascending=False)

    LOGGER.info(
        "Filtered {} gene sets down to {} after cutoffs"
        .format(len(vals), len(filtered_vals))
    )

    rows = max([len(filtered_vals) // cols, 1])
    scale = 6

    f, axes = plt.subplots(
        rows, cols,
        figsize=(scale * cols, scale * rows),
        squeeze=False,
        sharex=True,
        sharey=True,
    )
    axes = [i for j in axes for i in j]

    for (index, ax), (set_id, row) in zip(
        enumerate(axes),
        filtered_vals.iterrows(),
    ):
        ax.plot(row["cumscore"])

        for ind, hit in enumerate(row["hits"]):
            if hit:
                ax.axvline(ind, linestyle=":", alpha=.25)

        name = gene_sets.loc[set_id]["name"]
        name = name if len(name) < 40 else name[:40] + "..."

        ax.set_title(
            "{}\nhits: {} ES={:.2f}"
            .format(
                name,
                row["n_hits"],
                row["ES(S)"],
            ) + (
                ", p={:.2f}".format(
                    row["p-value"],
                ) if "p-value" in row.index else ""
            ) + (
                ", q={:.2f}".format(
                    row["q-value"],
                ) if "q-value" in row.index else ""
            ),
            fontsize=20,
        )
        ax.axhline(0, color="k")

        if index >= len(axes) - cols:
            ax.set_xlabel("Gene List Rank", fontsize=20)

        if index % cols == 0:
            ax.set_ylabel("ES(S)", fontsize=20)

        ax.set_ylim(-1, 1)

    return vals
