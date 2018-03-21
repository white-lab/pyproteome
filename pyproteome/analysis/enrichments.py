
from __future__ import division

from collections import defaultdict
from functools import partial
import logging
import multiprocessing

from adjustText.adjustText import adjust_text
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.utils import shuffle


LOGGER = logging.getLogger("pyproteome.enrichments")
MIN_PERIODS = 5
CORRELATION_METRICS = [
    "spearman",
    "pearson",
    "kendall",
    "fold",
    "zscore",
]


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
    ds, gene_sets,
    phenotype=None,
    p=1,
    metric="spearman",
    p_iter=1000,
):
    assert metric in CORRELATION_METRICS

    LOGGER.info(
        "Calculating ES(S, pi) for {} gene sets".format(len(gene_sets))
    )
    n_cpus = 6

    if metric in ["spearman", "pearson", "kendall"]:
        n_cpus = 4

    pool = multiprocessing.Pool(
        processes=n_cpus,
    )

    ess_dist = defaultdict(list)

    for ind, ess in enumerate(
        pool.imap_unordered(
            partial(
                _calc_essdist,
                ds=ds,
                gene_sets=gene_sets,
                p=p,
                metric=metric,
            ),
            [phenotype for _ in range(p_iter)],
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
    if nes > 0:
        return (
            nes_pi_pdf.sf(nes)
        ) / (
            nes_pdf.pdf(nes) + nes_pdf.sf(nes)
        )
    else:
        return (
            nes_pi_pdf.cdf(nes)
        ) / (
            nes_pdf.pdf(nes) + nes_pdf.cdf(nes)
        )


class PrPDF(object):
    def __init__(self, data):
        self.data = np.sort(data)

    def pdf(self, x):
        return 1 / self.data.shape[0]

    def cdf(self, x):
        return np.searchsorted(self.data, x, side="left") / self.data.shape[0]

    def sf(self, x):
        return 1 - (
            np.searchsorted(self.data, x, side="right") / self.data.shape[0]
        )


def estimate_pq(vals):
    vals = vals.copy()

    mask = vals["ES(S)"] > 0

    pos_ess = vals["ES(S)"][mask]
    neg_ess = vals["ES(S)"][~mask]

    pos_pi = vals["ES(S, pi)"].apply(np.array).apply(lambda x: x[x > 0])
    neg_pi = vals["ES(S, pi)"].apply(np.array).apply(lambda x: x[x < 0])

    pos_mean = pos_pi.apply(np.mean)
    neg_mean = neg_pi.apply(np.mean)

    pos_nes = (pos_ess / pos_mean)
    neg_nes = -(neg_ess / neg_mean)

    assert (pos_nes.isnull() | neg_nes.isnull()).all()
    assert ((~pos_nes.isnull()) | (~neg_nes.isnull())).all()

    pos_pi_nes = pos_pi / pos_mean
    neg_pi_nes = -neg_pi / neg_mean

    vals["NES(S)"] = pos_nes.fillna(neg_nes)
    vals["pos NES(S, pi)"] = pos_pi_nes
    vals["neg NES(S, pi)"] = neg_pi_nes

    pos_mat = (
        np.concatenate(pos_pi_nes.as_matrix())
        if pos_pi_nes.shape[0] > 0 else
        np.array([])
    )
    neg_mat = (
        np.concatenate(neg_pi_nes.as_matrix())
        if neg_pi_nes.shape[0] > 0 else
        np.array([])
    )

    f, ax = plt.subplots()
    ax.hist(
        pos_mat,
        bins=50,
        color='k',
        alpha=.5,
    )
    ax.hist(
        neg_mat,
        bins=50,
        color='k',
        alpha=.5,
    )
    ax.hist(
        vals["NES(S)"].as_matrix()
        if vals["NES(S)"].shape[0] > 0 else
        [],
        bins=50,
        color='r',
        alpha=.5,
    )

    LOGGER.info("Calculated pos, neg distributions")

    pos_pdf = PrPDF(pos_nes.dropna().as_matrix())
    neg_pdf = PrPDF(neg_nes.dropna().as_matrix())

    pos_pi_pdf = PrPDF(pos_mat)
    neg_pi_pdf = PrPDF(neg_mat)

    LOGGER.info("Normalized NES distributions")

    if vals.shape[0] > 0:
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
                pos_pdf if x["NES(S)"] > 0 else neg_pdf,
                pos_pi_pdf if x["NES(S)"] > 0 else neg_pi_pdf,
            ),
            axis=1,
        )

    LOGGER.info("Calculated p, q values")
    return vals


def _frac_true(x):
    return sum(x) / max([len(x), 1])


def absmax(i):
    return max([abs(min(i)), abs(max(i))])


def get_gene_changes(ds):
    LOGGER.info("Getting gene correlations")

    gene_changes = ds.psms[["ID", "Correlation"]]
    gene_changes.is_copy = None

    gene_changes["Abs Corr"] = gene_changes["Correlation"].apply(abs)
    gene_changes = gene_changes.sort_values(by="Abs Corr", ascending=True)
    gene_changes = gene_changes.drop_duplicates(subset="ID", keep="last")

    gene_changes = gene_changes.sort_values(by="Correlation", ascending=False)
    gene_changes = gene_changes.set_index(keys="ID")

    return gene_changes


def _calc_essdist(phen, ds=None, gene_sets=None, p=1, metric="spearman"):
    assert metric in CORRELATION_METRICS

    if metric in ["spearman", "pearson", "kendall"]:
        phen = _shuffle(phen)
    else:
        ds = ds.copy()
        ds.psms["Fold Change"] = _shuffle(ds.psms["Fold Change"])

    vals = enrichment_scores(
        ds,
        gene_sets,
        phenotype=phen,
        metric=metric,
        p=p,
        recorrelate=True,
        pval=False,
    )

    return vals["ES(S)"]


def correlate_phenotype(ds, phenotype=None, metric="spearman"):
    assert metric in CORRELATION_METRICS

    ds = ds.copy()

    if metric in ["spearman", "pearson", "kendall"]:
        ds.psms["Correlation"] = ds.psms.apply(
            lambda row:
            phenotype.corr(
                pd.to_numeric(row[phenotype.index]),
                method=metric,
                min_periods=MIN_PERIODS,
            ),
            axis=1,
        )
    else:
        new = ds.psms["Fold Change"]
        new = new.apply(np.log2)

        if metric in ["zscore"]:
            new = (new - new.mean()) / new.std()

        ds.psms["Correlation"] = new

    return ds


def calculate_es_s(gene_changes, gene_set, p=1):
    n = len(gene_changes)
    gene_set = [
        gene
        for gene in gene_set
        if gene in gene_changes.index
    ]
    hits = gene_changes.index.isin(gene_set)
    hit_list = gene_changes[hits].index.tolist()

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
    cumsum = np.cumsum(scores)

    # ess = cumsum.max() - cumsum.min()
    # ess *= np.sign(max(cumsum, key=abs))
    ess = max(cumsum, key=abs)

    return hits, hit_list, cumsum, ess


def enrichment_scores(
    ds, gene_sets,
    phenotype=None,
    p=1,
    pval=True,
    recorrelate=False,
    p_iter=1000,
    metric="spearman",
):
    """
    Here is the function for calculating gene set enrichment scores.
    It calculates the association of each gene with a given phenotype
    and then generates the ES(S) scores for a given gene set.
    """
    assert metric in CORRELATION_METRICS

    if recorrelate:
        ds = correlate_phenotype(ds, phenotype=phenotype, metric=metric)

    gene_changes = get_gene_changes(ds)

    cols = ["name", "cumscore", "ES(S)", "hits", "hit_list", "n_hits"]
    vals = pd.DataFrame(
        columns=cols,
    )

    for set_id, row in gene_sets.iterrows():
        hits, hit_list, cumscore, ess = calculate_es_s(
            gene_changes,
            row["set"],
            p=p,
        )
        vals = vals.append(
            pd.Series(
                [
                    row["name"],
                    cumscore,
                    ess,
                    hits,
                    hit_list,
                    len(hit_list),
                ],
                name=set_id,
                index=cols,
            )
        )

    if pval:
        vals = simulate_es_s_pi(
            vals, ds, gene_sets,
            phenotype=phenotype,
            p=p,
            metric=metric,
            p_iter=p_iter,
        )

        vals = estimate_pq(vals)
        vals = vals.sort_values("NES(S)", ascending=False)
    else:
        vals = vals.sort_values("ES(S)", ascending=False)

    return vals


def plot_nes(vals, max_pval=.1, max_qval=1):
    f, ax = plt.subplots()
    v = vals.copy()
    v = v.sort_values("NES(S)")
    mask = (
        (v["p-value"] < max_pval) &
        (v["q-value"] < max_qval)
    )

    nes = v["NES(S)"]

    pos = nes[mask]
    neg = nes[~mask]

    ind = np.arange(v.shape[0])
    pos_ind = ind[mask]
    neg_ind = ind[~mask]

    ax.scatter(
        x=pos_ind,
        y=pos,
        color="r",
    )
    ax.scatter(
        x=neg_ind,
        y=neg,
        color="k",
    )
    ax.set_ylabel("Normalized Enrichment Score (NES)")
    ax.set_xlim(
        left=ind.min() - 5,
        right=ind.max() + 5,
    )
    ax.set_ylim(
        bottom=nes.min() - .5,
        top=nes.max() + .5,
    )

    texts = []
    labels = []

    for m, x, (_, row) in zip(mask, ind, v.iterrows()):
        y = row["NES(S)"]
        labels.append(
            (x, y)
        )
        texts.append(
            ax.text(
                x=x,
                y=y,
                s=row["name"] if m else "",
                horizontalalignment='right',
            )
        )

    adjust_text(
        x=[i[0] for i in labels],
        y=[i[1] for i in labels],
        texts=texts,
        ax=ax,
        lim=500,
        force_text=0.5,
        force_points=0.1,
        arrowprops=dict(arrowstyle="-", relpos=(0, 0), lw=1),
        only_move={
            "points": "y",
            "text": "xy",
        }
    )

    return f, ax


def plot_correlations(gene_changes):
    """
    Plot the ranked list of correlations.
    """
    LOGGER.info("Plotting gene correlations")
    f, ax = plt.subplots()

    ax.plot(
        gene_changes["Correlation"].sort_values(ascending=False).tolist(),
    )
    ax.axhline(0, color="k")

    ax.set_xlabel("Gene List Rank", fontsize=20)
    ax.set_ylabel("Correlation", fontsize=20)

    return f, ax


def filter_gene_sets(gene_sets, ds, min_hits=10):
    LOGGER.info("Filtering gene sets")

    total_sets = gene_sets
    all_genes = set(ds["ID"])

    gene_sets = gene_sets[
        gene_sets["set"].apply(
            lambda x:
            len(x) < len(all_genes) and
            len(all_genes.intersection(x)) >= min_hits
        )
    ]

    LOGGER.info(
        "Filtered {} gene sets down to {} with {} or more genes present"
        .format(total_sets.shape[0], gene_sets.shape[0], min_hits)
    )

    return gene_sets


def filter_vals(
    vals,
    min_abs_score=.3,
    max_pval=1,
    max_qval=1,
):
    filtered_vals = vals[
        vals.apply(
            lambda x:
            abs(x["ES(S)"]) >= min_abs_score and (
                (x["p-value"] <= max_pval)
                if "p-value" in x.index else
                True
            ) and (
                (x["q-value"] <= max_qval)
                if "q-value" in x.index else
                True
            ),
            axis=1,
        )
    ]
    LOGGER.info(
        "Filtered {} gene sets down to {} after cutoffs"
        .format(len(vals), len(filtered_vals))
    )

    return filtered_vals


def plot_enrichment(
    vals, gene_sets,
    cols=5,
):
    rows = max([int(np.ceil(len(vals) / cols)), 1])
    scale = 6

    f, axes = plt.subplots(
        rows, cols,
        figsize=(scale * cols, scale * rows),
        squeeze=False,
        sharex=True,
        sharey=True,
    )
    axes = [i for j in axes for i in j]

    nes = "NES(S)" if "NES(S)" in vals.columns else "ES(S)"

    for (index, ax), (set_id, row) in zip(
        enumerate(axes),
        vals.iterrows(),
    ):
        ax.plot(row["cumscore"])

        for ind, hit in enumerate(row["hits"]):
            if hit:
                ax.axvline(ind, linestyle=":", alpha=.25)

        name = gene_sets.loc[set_id]["name"]
        name = name if len(name) < 35 else name[:35] + "..."

        ax.set_title(
            "{}\nhits: {} {}={:.2f}"
            .format(
                name,
                row["n_hits"],
                nes.split("(")[0],
                row[nes],
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

    return f, axes


# Here's a helper function to generate the plots seen below
def plot_gsea(
    ds, gene_sets,
    phenotype=None,
    p=1,
    cols=5,
    min_hits=10,
    min_abs_score=.3,
    max_pval=1,
    max_qval=1,
    pval=False,
    p_iter=1000,
    metric="spearman",
):
    assert metric in CORRELATION_METRICS

    gene_changes = get_gene_changes(ds)

    figs = ()

    figs += plot_correlations(gene_changes)[0],

    gene_sets = filter_gene_sets(
        gene_sets, ds,
        min_hits=min_hits,
    )

    vals = enrichment_scores(
        ds,
        gene_sets,
        phenotype=phenotype,
        p=p,
        pval=pval,
        p_iter=p_iter,
        metric=metric,
    )

    if vals.shape[0] < 1:
        return vals, figs

    figs += plot_enrichment(
        filter_vals(
            vals,
            min_abs_score=min_abs_score,
            max_pval=max_pval,
            max_qval=max_qval,
        ),
        gene_sets,
        cols=cols,
    )[0],

    if pval:
        figs += plot_nes(
            vals,
            max_pval=max_pval,
            max_qval=max_qval,
        )[0],

    return vals, figs
