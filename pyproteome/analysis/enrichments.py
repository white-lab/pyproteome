
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

# Here is the function for calculating gene set enrichment scores.
# It calculates the association of each gene with a given phenotype
# and then generates the ES(S) scores for a given gene set.


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


def enrichment_scores(
    ds, gene_sets, phenotype,
    p=1, pval=True, recorrelate=False, p_iter=1000,
):
    if recorrelate:
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

    gene_changes = _get_changes(ds)

    cols = ["cumscore", "ES(S)", "hits", "n_hits"]
    vals = pd.DataFrame(
        columns=cols,
    )

    for set_id, row in gene_sets.iterrows():
        n = len(gene_changes)
        gene_set = [
            gene
            for gene in row["set"]
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
        cumscore = np.cumsum(scores)
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
        LOGGER.info(
            "Calculating p-value for {} gene sets".format(len(gene_sets))
        )

        pool = multiprocessing.Pool(
            processes=6,
        )

        p_df = pd.DataFrame(
            columns=["NES(S)", "NES(S, pi)"],
            index=vals.index.copy(),
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

        LOGGER.info("Calculated pos, neg means")

        neg = vals["ES(S, pi)"].apply(
            lambda x:
            np.mean([i for i in x if i < 0] or [1])
        )
        pos = vals["ES(S, pi)"].apply(
            lambda x:
            np.mean([i for i in x if i > 0] or [1])
        )

        LOGGER.info("Calculated pos, neg distributions")

        def _norm_data(ess, key):
            div = np.array([(neg if i < 0 else pos)[key] for i in ess])
            return ess / div

        # p_df["NES(S)"] = vals.apply(
        #     lambda x:
        #     x["ES(S)"] / (neg if x["ES(S)"] < 0 else pos)[x.name],
        #     axis=1,
        # )
        #
        # p_df["NES(S, pi)"] = p_df.apply(
        #     lambda x:
        #     _norm_data(x["ES(S, pi)"], x.name),
        #     axis=1,
        # )
        LOGGER.info("Normalized NES distributions")

        vals["p-value"] = vals.apply(
            lambda x:
            _frac_true(
                _calc_p(x["ES(S)"], x["ES(S, pi)"])
            ),
            axis=1,
        )

        # vals["q-value"] = p_df.apply(
        #     lambda x:
        #     _frac_true(
        #         max(x["NES(S)"]) < max(x["NES(S, pi)"])
        #     ),
        #     axis=1,
        # )

        LOGGER.info("Calculated p, q values")

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
                ax.axvline(ind, linestyle="--", alpha=.25)

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
