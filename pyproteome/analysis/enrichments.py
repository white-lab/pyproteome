
from collections import OrderedDict
import logging
import numpy as np
from matplotlib import pyplot as plt

LOGGER = logging.getLogger("pyproteome.enrichments")

# Here is the function for calculating gene set enrichment scores.
# It calculates the association of each gene with a given phenotype
# and then generates the ES(S) scores for a given gene set.


def enrichment_scores(gene_changes, gene_set, p=1):
    if len(gene_changes) < len(gene_set):
        return None

    ranked_genes = sorted(
        gene_changes.items(),
        key=lambda x: x[1],
        reverse=True,
    )

    n = len(gene_changes)
    n_h = sum(
        gene in gene_changes
        for gene in gene_set
    )
    n_r = sum(
        abs(gene_changes.get(gene, 0)) ** p
        for gene in gene_set
    )
    hits = [
        gene in gene_set
        for gene, _ in ranked_genes
    ]

    scores = [0] + [
        ((abs(val) ** p) / n_r)
        if hit else
        (-1 / (n - n_h))
        for hit, (_, val) in zip(hits, ranked_genes)
    ]

    return np.cumsum(scores), hits


# Here's a helper function to generate the plots seen below
def plot_enrichment(
    sorted_set, gene_sets,
    p=1,
    cols=5,
    min_hits=10,
    min_abs_score=0.5,
):
    # Plot the ranked list of correlations
    f, ax = plt.subplots()
    ax.plot(sorted(
        sorted_set.values(),
        reverse=True,
    ))
    ax.axhline(0, color="k")
    ax.set_xlabel("Gene List Rank", fontsize=20)
    ax.set_ylabel("Correlation", fontsize=20)

    gene_sets = gene_sets[
        gene_sets["set"].apply(lambda x: len(x) < len(sorted_set))
    ]

    vals = OrderedDict()

    for index, row in gene_sets.iterrows():
        ess, hits = enrichment_scores(sorted_set, row["set"], p=p)
        vals[index] = (ess, hits)

    filtered_vals = sorted(
        [
            (index, (ess, hits))
            for index, (ess, hits) in vals.items()
            if sum(hits) >= min_hits and
            max([abs(max(ess)), abs(min(ess))]) >= min_abs_score
        ],
        key=lambda x: max([abs(max(x[1][0])), abs(min(x[1][0]))]),
        reverse=True,
    )

    LOGGER.info(
        "Filtered {} gene sets down to {} after cutoffs"
        .format(len(vals), len(filtered_vals))
    )

    rows = len(filtered_vals) // cols
    scale = 6

    f, axes = plt.subplots(
        rows, cols,
        figsize=(scale * cols, scale * rows),
        squeeze=False,
        sharex=True,
        sharey=True,
    )
    axes = [i for j in axes for i in j]

    for (index, ax), (set_id, (ess, hits)) in zip(
        enumerate(axes),
        filtered_vals,
    ):
        ax.plot(ess)

        for ind, hit in enumerate(hits):
            if hit:
                ax.axvline(ind, linestyle="--", alpha=.25)

        name = gene_sets.loc[set_id]["name"]
        name = name if len(name) < 40 else name[:40] + "..."

        ax.set_title(
            "{}\nhits: {} min/max ES: {:.2f}, {:.2f}"
            .format(name, sum(hits), min(ess), max(ess)),
            fontsize=20,
        )
        ax.axhline(0, color="k")

        if index >= len(axes) - cols:
            ax.set_xlabel("Gene List Rank", fontsize=20)

        if index % cols == 0:
            ax.set_ylabel("ES(S)", fontsize=20)

        ax.set_ylim(-1, 1)

    return vals


def scores(sorted_set, gene_set, permute_n=1000, p=1):
    np.random.seed(0)
    es_s = max(
        enrichment_scores(
            sorted_set,
            gene_set,
            p=p,
        ),
        key=lambda x: abs(x),
    )
    es_s_pi = np.array([
        max(
            enrichment_scores(
                dict(zip(
                    np.random.permutation(list(sorted_set.keys())),
                    sorted_set.values()
                )),
                gene_set,
                p=p,
            ),
            key=lambda x: abs(x),
        )
        for i in range(permute_n)
    ])

    # Only select the subset of the binomial distribution that we care
    # about. Not 100% sure if this is kosher, statistically.
    if es_s >= 0:
        es_s_pi = es_s_pi[es_s_pi >= 0]
    else:
        es_s_pi = es_s_pi[es_s_pi <= 0]

    nes_s = es_s / np.mean(es_s_pi)
    nes_s_pi = es_s_pi / np.mean(es_s_pi)

    return nes_s, nes_s_pi


def plot_permutations(sorted_set, gene_sets, cols=1):
    rows = len(gene_sets) // cols

    f, axes = plt.subplots(
        rows, cols,
        figsize=(2 * cols, 2 * rows),
        squeeze=False,
        sharey=True,
    )
    f.set_size_inches(cols * 3.5, 4)
    axes = [i for j in axes for i in j]

    # Compute the distributions of all (S, \pi) and (S) for
    # later use in calculating the FDR q-value
    fdr_nes_s, fdr_nes_s_pi = np.array([]), np.array([])

    for gene_set in gene_sets:
        nes_s, nes_s_pi = scores(sorted_set, gene_set)
        fdr_nes_s = np.append(fdr_nes_s, [nes_s])
        fdr_nes_s_pi = np.append(fdr_nes_s_pi, nes_s_pi)

    # Plot each distribution and compute the true p and q-values
    for index, ax, gene_set in zip(range(len(axes)), axes, gene_sets):
        nes_s, nes_s_pi = scores(sorted_set, gene_set)
        p_value = sum(nes_s_pi >= nes_s) / nes_s_pi.size

        q_nominator = sum(fdr_nes_s_pi >= nes_s) / fdr_nes_s_pi.size
        q_denomator = sum(fdr_nes_s >= nes_s) / fdr_nes_s.size

        q_value = q_nominator / q_denomator

        ax.hist(nes_s_pi, normed=1)
        ax.axvline(nes_s, color="red")
        ax.set_title(
            (
                "Permuted phenotype labels\nGene set: {}\n"
                "p-value={:.3f}, FDR q-value={:.3f}"
            ).format("gene_set", p_value, q_value)
        )

        if index >= len(axes) - cols:
            ax.set_xlabel(r"NES(S, $\pi$)", fontsize=20)

        if index % cols == 0:
            ax.set_ylabel("ES(S)", fontsize=20)

        ax.set_ylim(-1, 1)
