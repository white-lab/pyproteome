
from collections import Counter
import logging
import os
import re

from matplotlib import transforms
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from matplotlib.text import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties
import numpy as np
from scipy import stats

import pyproteome as pyp
from . import motif, plogo


BASES = list('ACDEFGHIKLMNPQRSTVWY')
GLOBSCALE = 1.4
LETTERS = {
    base: TextPath(
        (-0.303, 0),
        base,
        size=1,
        prop=FontProperties(family='monospace', weight='bold'),
    )
    for base in BASES
}
LETTERS['Q'] = TextPath(
    (-0.303, .11),
    'Q',
    size=1,
    prop=FontProperties(family='monospace', weight='bold'),
)
LETTERS['G'] = TextPath(
    (-0.303, .01),
    'G',
    size=1,
    prop=FontProperties(family='monospace', weight='bold'),
)
LETTER_YSCALE = {
    'Q': .84,
    'G': .95,
}

COLORS_SCHEME = {
    i: 'black'
    for i in BASES
}
COLORS_SCHEME.update({
    'C': '#BEB86B',
    'D': '#800000',
    'E': '#800000',
    'F': '#6F6F6F',
    'G': '#155939',
    'H': '#142B4F',
    'K': '#142B4F',
    'R': '#142B4F',
    'N': '#A97C50',
    'P': '#1C5E3F',
    'Q': '#A97C50',
    'S': '#4A79A5',
    'T': '#4A79A5',
    'L': '#000000',
    'A': '#000000',
    'I': '#000000',
    'M': '#000000',
    'V': '#000000',
    'W': '#000000',
    'Y': '#6F6F6F',
})
LOGGER = logging.getLogger('pyp.motifs.logo')


def _letterAt(letter, x, y, alpha=1, xscale=1, yscale=1, ax=None):
    text = LETTERS[letter]

    yscale *= LETTER_YSCALE.get(letter, .98)

    t = transforms.Affine2D().scale(
        xscale * GLOBSCALE, yscale * GLOBSCALE
    ) + transforms.Affine2D().translate(x, y) + ax.transData

    p = PathPatch(
        text,
        lw=0,
        fc=COLORS_SCHEME[letter],
        alpha=alpha,
        transform=t,
    )

    if ax is not None:
        ax.add_artist(p)

    return p


def _calc_score(
    fore_hit_size, fore_size, back_hit_size, back_size,
    prob_fn=None,
):
    if prob_fn is None:
        prob_fn = 'hypergeom'

    assert prob_fn in ['hypergeom', 'binom']

    if back_hit_size <= 0:
        return 0

    k = fore_hit_size
    n = fore_size
    K = back_hit_size
    N = back_size
    p = K / N

    if prob_fn == 'hypergeom':
        binomial = stats.hypergeom(N, K, n)
    else:
        binomial = stats.binom(n, p)

    pr_gt_k = binomial.sf(k - 1)
    pr_lt_k = binomial.cdf(k)

    if pr_lt_k <= 0:
        return -200
    elif pr_gt_k <= 0:
        return 200
    else:
        return -np.log10(pr_gt_k / pr_lt_k)


def _calc_scores(bases, fore, back, p=0.05, prob_fn=None):
    length = len(back[0])
    fore_counts = [
        Counter(i[pos] for i in fore)
        for pos in range(length)
    ]
    back_counts = [
        Counter(i[pos] for i in back)
        for pos in range(length)
    ]
    return {
        base: [
            _calc_score(
                fore_counts[pos][base],
                len(fore),
                back_counts[pos][base],
                len(back),
                prob_fn=prob_fn,
            )
            for pos in range(length)
        ]
        for base in bases
    }, _calc_hline(back_counts, p=p)


def _calc_hline(back_counts, p=0.05):
    '''
    Calculate the significance cutoff using multiple-hypothesis correction.

    Parameters
    ----------
    back_counts : collections.Counter of str, int
        Frequency of residues found in the background set.
    p : float, optional

    Returns
    -------
    float
        Signficance cutoff in log-odds space.
    '''
    num_calc = sum(
        1
        for counts in back_counts
        for _, count in counts.items()
        if count > 0
    )
    alpha = p / num_calc
    return abs(np.log10(alpha / (1 - alpha)))


def make_logo(data, f, **kwargs):
    '''
    Create a logo from a pyproteome data set using a given filter to define
    the foreground set.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.DataSet`
    f : dict
        Filter passed to :func:`pyproteome.data_sets.DataSet.filter` to define the foreground set.
    kwargs
        Arguments passed on to :func:`.logo`

    Returns
    -------
    fig, axes
    '''
    LOGGER.info('Generating motif logo')

    nmer_args = motif.get_nmer_args(kwargs)

    fore = [
        n.upper()
        for n in motif.generate_n_mers(
            data.filter(f)['Sequence'],
            **nmer_args
        )
    ]

    back = [
        n.upper()
        for n in motif.generate_n_mers(
            data['Sequence'],
            **nmer_args
        )
    ]
    title = kwargs.pop('title', plogo.format_title(data=data, f=f))
    fig, ax = logo(
        fore, back,
        title=title,
        **kwargs
    )

    return fig, ax


def _draw_logo(
    scores,
    ax,
    p_line=None,
    title=None,
    ytitle='',
    width=10,
    height=6,
    fade_power=1,
    low_res_cutoff=0,
    show_title=True,
    show_ylabel=True,
    minmaxy=None,
):
    length = len(list(scores.values())[0])
    left_margin = (
        .15 / width * 5
    )
    if show_ylabel:
        left_margin += .02

    ax.add_patch(
        patches.Rectangle(
            (left_margin, 0.01),
            .998 - left_margin,
            .98,
            fill=False,
            linewidth=1,
            edgecolor='k',
            zorder=10,
        )
    )
    ax.add_patch(
        patches.Rectangle(
            (left_margin, .46),
            .9985 - left_margin,
            .08,
            fill=False,
            linewidth=1,
            edgecolor='k',
            zorder=10,
        )
    )
    # ax.add_patch(
    #     patches.Rectangle(
    #         (left_margin, .5),
    #         .9985 - left_margin,
    #         .001,
    #         fill=False,
    #         linewidth=1,
    #         edgecolor='k',
    #         zorder=10,
    #     )
    # )

    axes = (
        ax.inset_axes([
            left_margin, .54,
            1 - left_margin, .46,
        ]),
        ax.inset_axes([
            left_margin, 0,
            1 - left_margin, .46,
        ])
    )
    yax = ax.inset_axes([
        0, 0,
        1, 1,
    ])
    xwidth = (1 - left_margin) / length
    xpad = xwidth / 2
    xax = ax.inset_axes([
        left_margin + xpad, 0.52,
        xwidth * (length - 1), .11,
    ])
    yax.patch.set_alpha(0)
    xax.patch.set_alpha(0)

    if p_line is not None:
        axes[0].axhline(p_line, color='red')
        axes[1].axhline(-p_line, color='red')

        miny, maxy = -p_line, p_line
    else:
        miny, maxy = 0, 0

    x = 1

    yax.xaxis.set_ticks([])
    yax.yaxis.set_ticks([])
    xax.yaxis.set_ticks([])
    xax.spines['bottom'].set_position(('data', 0))
    xax.set_ylim(bottom=-2, top=2.4)

    for ax in (yax, xax) + axes:
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.spines['right'].set_color('none')

    if show_title:
        yax.set_title(title)

    xax.set_xticks(
        range(0, length),
    )
    y_offset = (
        76 * np.power(xax.get_window_extent().height, -1.453)
    ) - .4
    y_offset = -.15

    xax.set_xticklabels(
        [
            '{:+d}'.format(i) if i != 0 else '0'
            for i in range(-(length - 1) // 2, (length - 1) // 2 + 1)
        ],
        va='center',
        ha='center',
        y=y_offset,
        fontsize=8,
    )

    for i in range(0, length):
        base_scores = [(b, scores[b][i]) for b in BASES]
        base_scores = (
            sorted([i for i in base_scores if i[1] < 0], key=lambda t: -t[1]) +
            sorted([i for i in base_scores if i[1] >= 0], key=lambda t: -t[1])
        )
        base_scores = [
            i
            for i in base_scores
            if abs(i[1]) >= (p_line or 0) * low_res_cutoff
        ]

        y = sum(i[1] for i in base_scores if i[1] < 0)
        miny = min(miny, y)

        for base, score in base_scores:
            _letterAt(
                base, x, y,
                alpha=min([1, abs(score / (p_line or 1))]) ** fade_power,
                xscale=1.2,
                yscale=abs(score),
                ax=axes[1 if score < 0 else 0],
            )
            y += abs(score)

        x += 1
        maxy = max(maxy, y)

    if minmaxy is None:
        minmaxy = max(abs(i) for i in [miny, maxy])

    for ind, ax in enumerate(axes):
        ax.set_xlim(
            left=.5,
            right=x - .5,
        )
        ax.set_ylim(
            bottom=-1.05 * minmaxy if ind == 1 else 0,
            top=1.05 * minmaxy if ind == 0 else 0,
        )
        ax.set_xticks([])

        spacing = minmaxy // 3
        if spacing != 0:
            ax.set_yticks(
                [
                    i
                    for i in np.arange(
                        spacing if ind == 0 else -spacing,
                        (spacing + 1) * (3 if ind == 0 else -3),
                        spacing * (1 if ind == 0 else -1)
                    )
                    if abs(i) >= abs(p_line)
                ],
            )
        else:
            ax.set_yticks(
                np.arange(
                    0,
                    minmaxy if ind == 0 else -minmaxy,
                    1 if ind == 0 else -1,
                )
            )

        ax.set_yticklabels(
            ax.get_yticks(),
        )
        if show_ylabel:
            yax.set_ylabel(
                ytitle,
            )

    return (yax, xax,) + axes, minmaxy


def logo(
    fore, back,
    ax=None,
    title='',
    width=12,
    height=8,
    p=0.05,
    fade_power=1,
    low_res_cutoff=0,
    prob_fn=None,
    show_title=True,
    show_ylabel=True,
    show_n=True,
    minmaxy=None,
):
    '''
    Generate a sequence logo locally using pLogo's enrichment score.

    Parameters
    ----------
    fore : list of str
    back : list of str
    title : str, optional
    p : float, optional
        p-value to use for residue significance cutoff. This value is corrected
        for multiple-hypothesis testing before being used.
    fade_power : float, optional
        Set transparency of residues with scores below p to:
        (score / p) ** fade_power.
    low_res_cutoff : float, optional
        Hide residues with scores below p * low_res_cutoff.
    prob_fn : str, optional
        Probability function to use for calculating enrichment. Either
        'hypergeom' or 'binom'. The default, hypergeom, is more accurate but
        more computationally expensive.

    Returns
    -------
    fig, axes
    '''
    if len(back) == 0:
        return None, None

    length = len(back[0])
    assert length > 0
    assert (
        all(len(i) == len(back[0]) for i in fore) and
        all(len(i) == len(back[0]) for i in back)
    )

    rel_info, p_line = _calc_scores(
        BASES, fore, back,
        p=p,
        prob_fn=prob_fn,
    )

    if ax is None:
        _, ax = plt.subplots(figsize=(width / 2, height / 2))
        ax.axis('off')

    axes, minmaxy = _draw_logo(
        scores=rel_info,
        p_line=p_line,
        title=title,
        ytitle='log odds',
        width=width,
        height=height,
        fade_power=fade_power,
        low_res_cutoff=low_res_cutoff,
        show_title=show_title,
        show_ylabel=show_ylabel,
        ax=ax,
        minmaxy=minmaxy,
    )

    if show_n:
        axes[3].text(
            length + .4,
            -minmaxy,
            'n(fg) = {}\nn(bg) = {}'.format(len(fore), len(back)),
            color='darkred',
            fontsize=18,
            ha='right',
            va='bottom',
        )

    return ax.get_figure(), axes
