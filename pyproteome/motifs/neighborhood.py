
from matplotlib import pyplot as plt
from scipy import stats

from . import motif, plogo


def enriched_neighborhood(
    data,
    f,
    residues,
    nmer_length=7,
    count_cutoff=2,
    mods=None,
):
    '''
    Calculates the hypergeometric enrichment value for the number of
    adjacent residues within a given window around all modification sites
    in a data set.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.data_set.DataSet`
    f : dict or list of dict
    residues : list of str
    nmer_length : int, optional
    count_cutoff : int, optional
    mods : str or list of str

    Returns
    -------
    f : :class:`matplotlib.figure.Figure`
    ax : :class:`matplotlib.axes.Axes`
    pval : float
        P-value, calculated with :class:`scipy.stats.hypergeom`.
    K : int
        Number of sequences with # residues > count_cutoff in background list.
    N : int
        Size of the background list of sequences.
    k : int
        Number of sequences with # residues > count_cutoff in foreground list.
    n : int
        Size of the foreground list of sequences.
    '''
    if mods is None:
        mods = [(None, 'Phospho')]
    background = motif.generate_n_mers(
        data['Sequence'],
        mods=mods,
        n=nmer_length,
        all_matches=False,
    )
    foreground = motif.generate_n_mers(
        data.filter(f)['Sequence'],
        mods=mods,
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
            color='green',
            bins=range(0, nmer_length, 1),
            label='background',
        )

    if foreground:
        ax.hist(
            [
                sum(i.count(j) for j in residues)
                for i in foreground
            ],
            density=True,
            alpha=0.7,
            color='orange',
            bins=range(0, nmer_length, 1),
            label=plogo.format_title(f=f),
        )
        ax.legend()

    ax.set_ylabel('Frequency')

    return fig, ax, pval, K, N, k, n
