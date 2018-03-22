
from matplotlib import pyplot as plt
from scipy import stats

from . import motif, plogo


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
        data.filter(f)["Sequence"],
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

    fig, ax = plt.subplots(figsize=(4, 4))

    if background:
        ax.hist(
            [
                sum(i.count(j) for j in residues)
                for i in background
            ],
            density=True,
            alpha=0.5,
            color="green",
            bins=range(0, nmer_length, 1),
            label="background",
        )

    if foreground:
        ax.hist(
            [
                sum(i.count(j) for j in residues)
                for i in foreground
            ],
            density=True,
            alpha=0.7,
            color="orange",
            bins=range(0, nmer_length, 1),
            label=plogo.format_title(f=f),
        )
        ax.legend()

    ax.set_title(
        "SFK Substrate Enrichment\np = {:.2e}\nK={}, N={}, k={}, n={}"
        .format(pval, K, N, k, n)
    )
    ax.set_xlabel(
        "# of acidic residues within {} residues of pY"
        .format(nmer_length // 2)
    )
    ax.set_ylabel("Frequency")

    return fig, ax, pval
