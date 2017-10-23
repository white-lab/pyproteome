
from matplotlib import pyplot as plt
from scipy import stats

from . import motif


def enriched_neighborhood(
    data, f, residues,
    nmer_length=7, count_cutoff=2,
):
    background = motif.generate_n_mers(
        data["Sequence"],
        letter_mod_types=[(None, "Phospho")],
        n=nmer_length,
        all_matches=False,
    )
    foreground = motif.generate_n_mers(
        data.filter(**f)["Sequence"],
        letter_mod_types=[(None, "Phospho")],
        n=nmer_length,
        all_matches=False,
    )

    N = len(background)
    K = len([
        i
        for i in background
        if sum(i.count(j) for j in residues) >= count_cutoff
    ])
    n = len(foreground)
    k = len([
        i
        for i in foreground
        if sum(i.count(j) for j in residues) >= count_cutoff
    ])
    pval = stats.hypergeom(
        N,
        K,
        n,
    ).sf(
        min([k, n]) - 1
    )

    print("K={}, N={}, k={}, n={}".format(K, N, k, n))
    print("p-value={:.3e}".format(pval))

    fig, ax = plt.subplots(figsize=(4, 4))

    ax.hist(
        [
            sum(i.count(j) for j in residues)
            for i in motif.generate_n_mers(
                data["Sequence"],
                letter_mod_types=[(None, "Phospho")],
                n=nmer_length,
            )
        ],
        normed=True, alpha=0.5, color="green",
        bins=range(0, nmer_length, 1),
        label="background",
    )

    ax.hist(
        [
            sum(i.count(j) for j in residues)
            for i in motif.generate_n_mers(
                data.filter(**f)["Sequence"],
                letter_mod_types=[(None, "Phospho")],
                n=nmer_length,
                all_matches=False,
            )
        ],
        normed=True, alpha=0.7, color="orange",
        bins=range(0, nmer_length, 1),
        label="foreground",
    )

    ax.set_title("SFK Substrate Enrichment")
    ax.set_xlabel(
        "# of acidic residues within {} residues of pY"
        .format(nmer_length // 2)
    )
    ax.set_ylabel("Freuency")
    ax.legend()

    return f, ax, pval
