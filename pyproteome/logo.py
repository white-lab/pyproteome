

from collections import Counter
import scipy.stats as ss
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.text import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties

from . import motif, plogo


bases = list("ACDEFGHIKLMNPQRSTVWY")
fp = FontProperties(family="monospace", weight="bold")
GLOBSCALE = 1.4
LETTERS = {
    base: TextPath((-0.303, 0), base, size=1, prop=fp)
    for base in bases
}

COLORS_SCHEME = {
    i: "black"
    for i in bases
}
COLORS_SCHEME.update({
    "A": "#000000",
    "C": "#BEB86B",
    "D": "#800000",
    "E": "#800000",
    "F": "#6F6F6F",
    "G": "#155939",
    "H": "#142B4F",
    "I": "#000000",
    "K": "#142B4F",
    "L": "#000000",
    "M": "#000000",
    "N": "#A97C50",
    "P": "#1C5E3F",
    "Q": "#A97C50",
    "R": "#142B4F",
    "S": "#4A79A5",
    "T": "#4A79A5",
    "V": "#000000",
    "W": "#000000",
    "Y": "#6F6F6F",
})


def letterAt(letter, x, y, alpha=1, xscale=1, yscale=1, ax=None):
    text = LETTERS[letter]

    t = mpl.transforms.Affine2D().scale(
        xscale * GLOBSCALE, yscale * GLOBSCALE
    ) + mpl.transforms.Affine2D().translate(x, y) + ax.transData

    p = PathPatch(
        text,
        lw=0,
        fc=COLORS_SCHEME[letter],
        alpha=alpha,
        transform=t
    )

    if ax is not None:
        ax.add_artist(p)

    return p


def _calc_score(fore_hit_size, fore_size, back_hit_size, back_size, base, pos):
    if back_hit_size <= 0:
        return 0

    p = back_hit_size / back_size
    K = fore_hit_size
    N = fore_size

    binomial = ss.binom(N, p)

    pr_gt_k = binomial.sf(K)
    pr_lt_k = binomial.cdf(K)

    if pr_lt_k <= 0:
        return -200
    elif pr_gt_k <= 0:
        return 200
    else:
        return -np.log(pr_gt_k / pr_lt_k)


def _calc_scores(bases, fore, back):
    length = len(back[0])
    # print(fore, len(fore))
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
                base, pos,
            )
            for pos in range(length)
        ]
        for base in bases
    }, _calc_hline(back_counts)


def _calc_hline(back_counts):
    num_calc = sum(
        1
        for counts in back_counts
        for _, count in counts.items()
        if count > 0
    )
    alpha = 0.05 / num_calc
    return abs(np.log(alpha / (1 - alpha)))


def make_logo(data, f, **kwargs):
    nmer_args = motif.get_nmer_args(kwargs)
    low_res_cutoff = kwargs.get("low_res_cutoff", None)

    fore = [
        n.upper()
        for n in motif.generate_n_mers(
            data.filter(**f)["Sequence"],
            **nmer_args
        )
    ]
    back = [
        n.upper()
        for n in motif.generate_n_mers(
            data["Sequence"],
            **nmer_args
        )
    ]

    length = len(back[0])
    assert length > 0

    fig, ax = plt.subplots(figsize=(12, 8))

    rel_info, p_line = _calc_scores(bases, fore, back)

    ax.axhline(0, color="black")
    ax.axhline(p_line, color="red")
    ax.axhline(-p_line, color="red")

    miny, maxy = -p_line, p_line
    x = 1

    for i in range(0, length):
        scores = [(b, rel_info[b][i]) for b in bases]
        scores = (
            sorted([i for i in scores if i[1] < 0], key=lambda t: -t[1]) +
            sorted([i for i in scores if i[1] >= 0], key=lambda t: -t[1])
        )
        if low_res_cutoff:
            scores = [
                i
                for i in scores
                if abs(i[1]) >= p_line * low_res_cutoff
            ]

        y = sum(i[1] for i in scores if i[1] < 0)
        miny = min(miny, y)

        for base, score in scores:
            letterAt(
                base, x, y,
                alpha=min([1, abs(score / p_line)]),
                xscale=1.2,
                yscale=abs(score),
                ax=ax,
            )
            y += abs(score)

        x += 1
        maxy = max(maxy, y)

    ax.set_xlim((0, x))
    ax.set_ylim((miny * 1.05, maxy * 1.05))

    ax.set_xticks(range(1, x))
    ax.set_yticklabels(ax.get_yticks(), fontsize=16)
    ax.set_xticklabels([
        "{:+d}".format(i) if i != 0 else "0"
        for i in range(-int(np.floor(length / 2)), int(np.ceil(length / 2)))
    ], fontsize=16)
    ax.set_ylabel("log odds of the binomial probability", fontsize=20)

    ax.set_title(plogo.format_title(data, f), fontsize=32)

    return fig, ax
