# Here is the function for calculating gene set enrichment scores.
# It calculates the association of each gene with a given phenotype
# and then generates the ES(S) scores for a given gene set.

def enrichment_scores(gene_changes, gene_set, p=1):
    assert len(gene_changes) >= len(gene_set)

    ranked_genes = sorted(gene_changes.items(), key=lambda x: x[1], reverse=True)

    n = len(gene_changes)
    n_h = len(gene_set)
    n_r = sum(
        abs(gene_changes.get(gene, 0)) ** p 
        for gene in gene_set
    )

    scores = [0] + [
        (abs(val) ** p) / n_r
        if gene in gene_set
        else
        - 1 / (n - n_h)
        for gene, val in ranked_genes
    ]

    return np.cumsum(scores)

# Here's a helper function to generate the plots seen below
def plot_enrichment(sorted_set, gene_sets, p=1, cols=1):
    f, axes = plt.subplots(len(gene_sets) // cols, cols, squeeze=False, sharey=True)
    f.set_size_inches(cols * 3.5, 4)
    axes = [i for j in axes for i in j]
    
    for index, ax, gene_set in zip(range(len(axes)), axes, gene_sets):
        ess = enrichment_scores(sorted_set, gene_set, p=p)
        ax.plot(ess)
        ax.set_title("{}, p={}\n(max: {:.2f}, min: {:.2f})"
                     .format("gene_set", p, max(ess), min(ess)))
        ax.axhline(0)
        ax.set_xlabel("Gene List Rank")
        
        if index % cols == 0:
            ax.set_ylabel("ES(S)")

def scores(sorted_set, gene_set, permute_n=1000):
    seed(0)
    es_s = max(
        enrichment_scores(
            sorted_set, 
            gene_set,
            p=1,
        ),
        key=lambda x: abs(x),
    )
    es_s_pi = np.array([
        max(
            enrichment_scores(
                dict(zip(np.random.permutation(list(sorted_set.keys())), sorted_set.values())),
                gene_set,
                p=1,
            ),
            key = lambda x: abs(x),
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
    f, axes = plt.subplots(len(gene_sets) // cols, cols, squeeze=False, sharey=True)
    f.set_size_inches(cols * 3.5, 4)
    axes = [i for j in axes for i in j]

    # Compute the distributions of all (S, \pi) and (S) for
    # later use in calculating the FDR q-value
    fdr_nes_s, fdr_nes_s_pi = array([]), array([])

    for gene_set in gene_sets:
        nes_s, nes_s_pi = scores(sorted_set, gene_set)
        fdr_nes_s = append(fdr_nes_s, [nes_s])
        fdr_nes_s_pi = append(fdr_nes_s_pi, nes_s_pi)

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
        ax.set_xlabel("NES(S, $\pi$)")
        
        if index % cols == 0:
            ax.set_ylabel("Frequency")

def load_gene_set(path):
    with open(path) as f:
        return [
            line.split("#")[0].strip().upper()
            for line in f
        ]

