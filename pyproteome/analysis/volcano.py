
from __future__ import absolute_import, division

import logging
import os
import re

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

import pyproteome as pyp


LOGGER = logging.getLogger("pyproteome.volcano")
MAX_VOLCANO_LABELS = 500


def _remove_lesser_dups(labels, compress_sym=False):
    if labels.shape[0] < 1:
        return labels

    labels["xy"] = labels.apply(lambda x: abs(x["x"]) + x["y"], axis=1)
    labels = labels.sort_values("xy", ascending=False)

    if compress_sym:
        labels = labels.drop_duplicates(subset="Label")
    else:
        labels = pd.concat([
            labels[labels["x"] >= 0].drop_duplicates(subset="Label"),
            labels[labels["x"] < 0].drop_duplicates(subset="Label"),
        ])

    return labels


def _get_name(proteins):
    genes = sorted(proteins.genes)
    common = ""
    sep = " / "

    if len(genes) > 1:
        common = os.path.commonprefix(genes)

        if common:
            sep = "/"

    return common + sep.join(i[len(common):] for i in genes)


def plot_volcano_labels(
    data,
    ax,
    upper_fold=None,
    lower_fold=None,
    p=None,
    fold_and_p=True,
    sequence_labels=False,
    options=None,
    show_duplicates=False,
    compress_sym=True,
    adjust=True,
    mods=None,
):
    """
    Plot labels on a volcano plot.
    """
    options = options or {}

    highlight = options.get('highlight', {})
    hide = options.get('hide', {})
    show = options.get('show', {})
    edgecolors = options.get('edgecolors', {})
    rename = options.get('rename', {})

    xminmax = ax.get_xlim()
    yminmax = ax.get_ylim()

    labels = data.psms.copy()

    labels["x"] = labels["Fold Change"]
    labels["y"] = labels["p-value"]

    labels = labels[
        (labels["x"] >= xminmax[0]) &
        (labels["x"] <= xminmax[1]) &
        (labels["y"] >= yminmax[0]) &
        (labels["y"] <= yminmax[1])
    ]

    if sequence_labels:
        labels["Label"] = labels["Sequence"].apply(str)
    elif not mods:
        labels["Label"] = labels["Proteins"].apply(_get_name)
    else:
        labels["Label"] = labels.apply(
            lambda x:
            " / ".join(
                [
                    "{} {}".format(
                        rename.get(gene, gene),
                        x["Modifications"].get_mods(
                            mods,
                        ).__str__(prot_index=index),
                    )
                    for index, gene in enumerate(x["Proteins"].genes)
                ]
            )
            if len(list(x["Modifications"].get_mods(mods))) > 0 else
            " / ".join(x["Proteins"].genes),
            axis=1,
        )

    def _get_names(row):
        names = [
            rename.get(row["Label"], row["Label"]),
            row["Label"],
        ]
        if "Sequence" in row:
            names += [str(row["Sequence"])]
        if "Proteins" in row:
            names += [
                j
                for i in row["Proteins"].genes
                for j in [i, rename.get(i, i)]
            ]

        return names

    labels["Names"] = labels.apply(_get_names, axis=1)

    labels = labels[
        labels["Names"].apply(lambda x: not any([i in hide for i in x]))
    ]
    labels = labels[
        labels["Names"].apply(lambda x: any([i in show for i in x])) | (
            (
                (labels["y"] >= p) &
                (
                    (labels["x"] >= upper_fold) |
                    (labels["x"] <= lower_fold)
                )
            )
            if fold_and_p else
            (
                (labels["y"] >= p) |
                (
                    (labels["x"] >= upper_fold) |
                    (labels["x"] <= lower_fold)
                )
            )
        )
    ]
    labels["Highlight"] = labels["Names"].apply(
        lambda x:
        any([
            i in highlight
            for i in x
        ])
    )
    labels["Label"] = labels["Names"].apply(
        lambda x:
        x[0]
    )

    def _get_txt_color(row):
        colors = [
            edgecolors.get(i)
            for i in row["Names"][:2]
            if i in edgecolors
        ]
        edgecolor = colors[0] if colors else None

        if edgecolor is None:
            gene_colors = [
                edgecolors.get(gene, None)
                for gene in row["Proteins"].genes
            ]
            if len(set(gene_colors)) == 1:
                edgecolor = gene_colors[0]

        return edgecolor

    labels["EdgeColor"] = labels.apply(
        _get_txt_color,
        axis=1,
    ) if labels.shape[0] > 0 else []

    labels = labels[
        [
            "x",
            "y",
            "Label",
            "EdgeColor",
            "Highlight",
        ]
    ].copy()

    if not show_duplicates:
        labels = _remove_lesser_dups(labels, compress_sym=compress_sym)

    # Position the labels
    txt_lim = 100 if mods else 12

    texts = [
        ax.text(
            x=row["x"],
            y=row["y"],
            s=row["Label"][:txt_lim] + (
                "..." if len(row["Label"]) > txt_lim else ""
            ),
            zorder=10,
            fontsize=20 if row["Highlight"] else 12,
            horizontalalignment=(
                'left' if row["x"] > 0 else "right"
            ),
            bbox=dict(
                alpha=1,
                linewidth=0.1,
                pad=.2,
                facecolor=row["EdgeColor"] or (
                    "#DDDDDD"
                    if edgecolors else
                    ("#BFEE90" if row["x"] > 0 else "#FFC1C1")
                ),
                zorder=1,
                # edgecolor="black",
                boxstyle="round",
            ),
        )
        for _, row in labels.iterrows()
    ]

    LOGGER.info("Plotting volcano labels for {} peptides".format(len(texts)))

    if adjust:
        texts = texts[:MAX_VOLCANO_LABELS]
        pyp.utils.adjust_text(
            x=[i._x for i in texts],
            y=[i._y for i in texts],
            texts=texts,
            ax=ax,
            lim=100,
            force_text=0.2,
            force_points=0.01,
            arrowprops=dict(
                arrowstyle="->",
                relpos=(0, 0),
                lw=1,
                zorder=1,
                color="k",
            ),
            only_move={
                "points": "y",
                "text": "xy",
            }
        )

    return labels


def plot_volcano(
    data,
    group_a=None,
    group_b=None,
    p=0.05,
    fold=1.25,
    xminmax=None,
    yminmax=None,
    title=None,
    ax=None,
    filename=None,
    folder_name=None,
    show_xlabel=True,
    show_ylabel=True,
    log2_fold=True,
    log10_p=True,
    bonferoni=False,
    **kwargs
):
    """
    Display a volcano plot of data.

    This plot inclues the fold-changes and p-values associated with said
    changes.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.DataSet`
    group_a : str or list of str, optional
    group_b : str or list of str, optional
    p : float, optional
    fold : float, optional
    folder_name : str, optional
    title : str, optional
    figsize : tuple of float, float
    kwargs : dict
        Arguments passed to :func:`.plot_volcano_labels`

    Returns
    -------
    f : :class:`matplotlib.figure.Figure`
    ax : :class:`matplotlib.axes.Axes`
    folder_name : str
    filename : str
    """
    data = data.copy()

    (channels_a, channels_b), (label_a, label_b), _ = data.get_groups(
        group_a=group_a,
        group_b=group_b,
    )

    if group_a and group_b:
        data.update_group_changes(group_a=group_a, group_b=group_b)

    folder_name = pyp.utils.make_folder(
        data=data,
        folder_name=folder_name,
        sub="Volcano",
    )

    if not filename:
        filename = re.sub("[ ></?]", "_", data.name) + ".png"

    if log10_p:
        p = -np.log10(p)
        data.psms["p-value"] = data.psms["p-value"].apply(
            lambda x: -np.log10(x)
        )

    if log2_fold:
        fold = np.log2(fold)
        data.psms["Fold Change"] = data.psms["Fold Change"].apply(
            lambda x: np.log2(x)
        )

    upper_fold = fold
    lower_fold = -upper_fold

    with pd.option_context('mode.use_inf_as_null', True):
        data.psms = data.psms.dropna(
            subset=["p-value", "Fold Change"],
            how="any",
        )

    if bonferoni:
        p += np.log10(data.shape[0])

    # Draw the figure
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))

    ax.scatter(
        data["Fold Change"],
        data["p-value"],
        s=5,
        c="grey",
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if not np.isnan(p):
        ax.axhline(
            p,
            color="r", linestyle="dashed", linewidth=0.5,
        )

    if (abs(fold) if log2_fold else abs(fold - 1)) > 0.01:
        ax.axvline(
            upper_fold,
            color="r",
            linestyle="dashed",
            linewidth=0.5,
        )
        ax.axvline(
            lower_fold,
            color="r",
            linestyle="dashed",
            linewidth=0.5,
        )

    ax.set_xticks(
        list(
            sorted(
                tick
                for tick in ax.get_xticks()
                if tick < lower_fold - .25 or tick > upper_fold + .25
            ) + [lower_fold, upper_fold]
        )
    )
    ax.set_yticks(
        list(
            sorted(
                [
                    tick
                    for tick in ax.get_yticks()
                    if "{}".format(np.power(1/10, tick)).strip("0.")[:1] in
                    ["1", "5"] and
                    tick > p
                ] + [p]
            )
        )
    )

    if xminmax:
        ax.set_xlim(left=xminmax[0], right=xminmax[1])
    else:
        ax.set_xlim(
            left=np.floor(min(data["Fold Change"] + [0]) * 2) / 2,
            right=np.ceil(max(data["Fold Change"] + [0]) * 2) / 2,
        )

    if yminmax:
        ax.set_ylim(bottom=yminmax[0], top=yminmax[1])
    else:
        ax.set_ylim(bottom=-0.1)

    ax.set_xticklabels(
        [
            "{:.3}".format(np.exp2(i) if log2_fold else i)
            for i in ax.get_xticks()
        ],
    )
    ax.set_yticklabels(
        [
            (
                "{:.3}" if i > 5e-3 else "{:.0e}"
            ).format(np.power(1/10, i) if log10_p else i)
            for i in ax.get_yticks()
        ],
    )

    if show_xlabel:
        max_len = 25
        ax.set_xlabel(
            "{} (n={}) / {} (n={})".format(
                label_a[:max_len] + ("..." if len(label_a) > max_len else ""),
                len(channels_a),
                label_b[:max_len] + ("..." if len(label_b) > max_len else ""),
                len(channels_b),
            ),
        )

    if show_ylabel:
        ax.set_ylabel(
            "p-value",
        )

    if title:
        ax.set_title(
            title,
        )

    plot_volcano_labels(
        data=data,
        ax=ax,
        upper_fold=upper_fold,
        lower_fold=lower_fold,
        p=p,
        **kwargs
    )

    fig = ax.get_figure()

    if filename:
        fig.savefig(
            os.path.join(folder_name, filename),
            bbox_inches="tight",
            dpi=pyp.DEFAULT_DPI,
            transparent=True,
        )

    return fig, ax, folder_name, filename


def plot_volcano_filtered(data, f, **kwargs):
    """
    Display a volcano plot, showing only peptides that are included by a given
    filter.

    All extra arguments will be passed directly to `plot_volcano`.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.DataSet`
    f : dict or list of dict
        Filters passed to data.filter().

    Returns
    -------
    f : :class:`matplotlib.figure.Figure`
    ax : :class:`matplotlib.axes.Axes`
    """
    data = data.copy()
    data.update_group_changes(
        group_a=kwargs.get("group_a", None),
        group_b=kwargs.get("group_b", None),
    )

    d = data.filter(f)

    changes = []
    pvals = []

    for _, row in data.psms.iterrows():
        row_pval = -np.log10(row["p-value"])
        row_change = np.log2(row["Fold Change"])

        if (
            np.isnan(row_pval) or
            np.isnan(row_change) or
            np.isinf(row_pval) or
            np.isinf(row_change)
        ):
            continue

        pvals.append(row_pval)
        changes.append(row_change)

    xminmax = kwargs.pop(
        "xminmax",
        (
            np.floor(min(changes + [0]) * 2) / 2,
            np.ceil(max(changes + [0]) * 2) / 2,
        ),
    )
    yminmax = kwargs.pop(
        "yminmax",
        (-0.1, np.ceil(max(pvals + [1]))),
    )

    f, ax, folder_name, filename = plot_volcano(
        d,
        xminmax=xminmax,
        yminmax=yminmax,
        **kwargs
    )

    if changes and pvals:
        changes, pvals = zip(*[
            (x, y)
            for x, y in zip(changes, pvals)
            if x > xminmax[0] and x < xminmax[1] and
            y > yminmax[0] and y < yminmax[1]
        ])

    ax.scatter(
        changes,
        pvals,
        s=5,
        zorder=0,
        c="grey",
        # c="lightblue",
        # alpha=0.3,
    )

    if filename:
        f.savefig(
            os.path.join(folder_name, filename),
            bbox_inches="tight",
            dpi=pyp.DEFAULT_DPI,
            transparent=True,
        )

    return f, ax
