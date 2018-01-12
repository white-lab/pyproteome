
from __future__ import division

import logging
import os
import re

from matplotlib import pyplot as plt
import numpy as np

# Misc extras
from adjustText.adjustText import adjust_text

from . import utils


LOGGER = logging.getLogger("pyproteome.volcano")


def _remove_lesser_dups(labels, compress_sym=False):
    new_labels = []

    for index, (p, change, label, color, highlight) in enumerate(labels):
        for o_index, (o_p, o_change, o_label, _, _) in enumerate(labels):
            if index == o_index:
                continue
            if label != o_label:
                continue
            if not compress_sym and (change < 0) != (o_change < 0):
                continue

            mul = 1 if change >= 0 else -1
            o_mul = 1 if o_change >= 0 else -1

            if (
                p + mul * change < o_p + o_mul * o_change
            ) or (
                (
                    p + mul * change <= o_p + o_mul * o_change
                ) and
                index < o_index
            ):
                break
        else:
            new_labels.append((p, change, label, color, highlight))

    return new_labels


def _get_color(txt, x, y):
    return "#BFEE90" if x > 0 else "#FFC1C1"


def _get_folder(data):
    folder_name = os.path.join(data.name, "Volcano")
    utils.make_folder(folder_name)
    return folder_name


def volcano_plot(
    data,
    group_a=None,
    group_b=None,
    p_cutoff=0.05,
    fold_cutoff=1.25,
    fold_and_p=True,
    xminmax=None, yminmax=None,
    options=None,
    folder_name=None, title=None,
    figsize=(12, 10),
    adjust=True,
    compress_dups=True,
    compress_sym=False,
    show_xlabel=True,
    show_ylabel=True,
    show_title=True,
    full_site_labels=False,
    sequence_labels=False,
    alpha=1,
    bonferoni=False,
):
    """
    Display a volcano plot of data.

    This plot inclues the fold-changes and p-values associated with said
    changes.

    Parameters
    ----------
    data : :class:`DataSet<pyproteome.data_sets.DataSet>`
    group_a : str or list of str, optional
    group_b : str or list of str, optional
    pval_cutoff : float, optional
    fold_cutoff : float, optional
    options: dict of (str, list), optional
    folder_name : str, optional
    title : str, optional
    figsize : tuple of float, float
    adjust : bool, optional
    compress_dups : bool, optional
    full_site_labels : bool, optional
    """
    (channels_a, channels_b), (label_a, label_b), _ = data.get_groups(
        group_a=group_a,
        group_b=group_b,
    )

    if group_a and group_b:
        data.update_group_changes(group_a=group_a, group_b=group_b)

    if not folder_name:
        folder_name = _get_folder(data)

    options = options or {}

    highlight = options.get('highlight', {})
    hide = options.get('hide', {})
    show = options.get('show', {})
    edgecolors = options.get('edgecolors', {})
    rename = options.get('rename', {})

    log_p_cutoff = -np.log10(p_cutoff)

    if not title:
        title = "{} - {}".format(
            data.name,
            data.enrichment,
        )

    if title:
        file_name = re.sub("[ ></?]", "_", title) + "_Volcano.png"

        if folder_name:
            file_name = os.path.join(folder_name, file_name)

    upper_fold = np.log2(fold_cutoff)
    lower_fold = -upper_fold

    rows = [
        row
        for _, row in data.psms.dropna(
            subset=["p-value", "Fold Change"],
        ).iterrows()
        if not (
            np.isnan(-np.log10(row["p-value"])) or
            np.isnan(np.log2(row["Fold Change"])) or
            np.isinf(-np.log10(row["p-value"])) or
            np.isinf(np.log2(row["Fold Change"]))
        )
    ]

    if bonferoni:
        p_cutoff = p_cutoff / len(rows)
        log_p_cutoff = -np.log10(p_cutoff)

    # Calculate the Fold-Change / p-values
    pvals = []
    changes = []
    colors = []
    labels = []

    for row in rows:
        color = "grey"
        row_pval = -np.log10(row["p-value"])
        row_change = np.log2(row["Fold Change"])

        if xminmax and (row_change < xminmax[0] or row_change > xminmax[1]):
            continue

        if sequence_labels:
            row_label = row["Sequence"].__str__()
        else:
            row_label = " / ".join(sorted(row["Proteins"].genes))
        old_row_label, old_re_row_label = None, None

        if (
            full_site_labels and
            len(list(row["Modifications"].skip_labels_iter())) > 0
        ):
            old_row_label = row_label
            old_re_row_label = rename.get(old_row_label, old_row_label)
            row_label = " / ".join(
                sorted(
                    "{} {}".format(
                        gene,
                        row["Modifications"].__str__(prot_index=index),
                    )
                    for index, gene in enumerate(row["Proteins"].genes)
                )
            )

        re_row_label = rename.get(row_label, row_label)
        edgecolor = edgecolors.get(
            re_row_label, edgecolors.get(
                row_label,
                edgecolors.get(old_row_label, None),
            )
        )

        if edgecolor is None:
            gene_colors = [
                edgecolors.get(gene, None)
                for gene in row["Proteins"].genes
            ]
            if len(set(gene_colors)) == 1:
                edgecolor = gene_colors[0]

        pvals.append(row_pval)
        changes.append(row_change)

        names = [
            row_label, re_row_label, old_row_label, old_re_row_label,
        ]

        if (
            any(
                i in show
                for i in names
            ) or ((
                row_pval > log_p_cutoff and
                (row_change > upper_fold or row_change < lower_fold)
            ) if fold_and_p else (
                row_pval > log_p_cutoff or
                (row_change > upper_fold or row_change < lower_fold)
            ))
        ):
            color = "blue"
            highlight_label = any(i in highlight for i in names)

            if (
                not any(i.strip() in hide for i in row_label.split("/")) and
                re_row_label not in hide
            ):
                labels.append((
                    row_pval, row_change,
                    re_row_label, edgecolor, highlight_label,
                ))

        colors.append(color)

    if compress_dups:
        labels = _remove_lesser_dups(labels, compress_sym=compress_sym)

    # Draw the figure
    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(changes, pvals, c=colors, alpha=alpha)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if xminmax:
        ax.set_xlim(xmin=xminmax[0], xmax=xminmax[1])
    else:
        ax.set_xlim(
            xmin=np.floor(min(changes) * 2) / 2,
            xmax=np.ceil(max(changes) * 2) / 2,
        )

    if yminmax:
        ax.set_ylim(bottom=yminmax[0], top=yminmax[1])
    else:
        ax.set_ylim(bottom=-0.1)

    ax.set_xticks(
        list(
            sorted(
                tick
                for tick in ax.get_xticks()
                if tick < lower_fold or tick > upper_fold
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
                    tick > log_p_cutoff
                ] + [log_p_cutoff]
            )
        )
    )
    ax.set_xticklabels(
        [
            "{:.3}".format(i)
            for i in np.exp2(ax.get_xticks())
        ],
        fontsize=20,
    )
    ax.set_yticklabels(
        [
            ("{:.3}" if i > 0.005 else "{:.0e}").format(i)
            for i in np.power(1/10, ax.get_yticks())
        ],
        fontsize=20,
    )

    if show_xlabel:
        ax.set_xlabel(
            "Fold Change {} (n={}) / {} (n={})".format(
                label_a,
                len(channels_a),
                label_b,
                len(channels_b),
            ),
            fontsize=20,
        )

    if show_ylabel:
        ax.set_ylabel(
            "p-value",
            fontsize=20,
        )

    if not np.isnan(log_p_cutoff):
        ax.axhline(
            log_p_cutoff,
            color="r", linestyle="dashed", linewidth=0.5,
        )

    if abs(fold_cutoff - 1) > 0.01:
        ax.axvline(upper_fold, color="r", linestyle="dashed", linewidth=0.5)
        ax.axvline(lower_fold, color="r", linestyle="dashed", linewidth=0.5)

    # Position the labels
    texts = []
    txt_lim = 100 if full_site_labels else 20

    for y, x, txt, edgecolor, highlight_label in labels:
        text = ax.text(
            x, y,
            txt[:txt_lim] + ("..." if len(txt) > txt_lim else ""),
        )

        if highlight_label:
            text.set_fontsize(20)

        text.set_bbox(
            dict(
                facecolor=_get_color(txt, x, y),
                alpha=1,
                linewidth=0.5 if not edgecolor else 3,
                edgecolor=edgecolor or "black",
                boxstyle="round",
            )
        )

        texts.append(text)

    if adjust:
        adjust_text(
            x=[i[0] for i in labels],
            y=[i[1] for i in labels],
            texts=texts,
            ax=ax,
            lim=500,
            force_text=0.5,
            force_points=0.01,
            arrowprops=dict(arrowstyle="->", relpos=(0, 0), lw=1),
            only_move={
                "points": "y",
                "text": "xy",
            }
        )

    if title:
        if show_title:
            ax.set_title(
                title,
                fontsize=20,
            )

        fig.savefig(
            file_name,
            bbox_inches="tight", dpi=300,
            transparent=True,
        )
        fig.savefig(
            os.path.splitext(file_name)[0] + ".svg",
            bbox_inches="tight", dpi=300,
            transparent=True,
        )

    return fig, ax


def plot_volcano_filtered(data, f, **kwargs):
    data = data.copy()
    data.update_group_changes(
        group_a=kwargs.get("group_a", None),
        group_b=kwargs.get("group_b", None),
    )

    d = data.filter(**f)

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
            np.floor(min(changes) * 2) / 2,
            np.ceil(max(changes) * 2) / 2,
        ),
    )
    yminmax = kwargs.pop(
        "yminmax",
        (-0.1, np.ceil(max(pvals))),
    )

    f, ax = volcano_plot(
        d,
        xminmax=xminmax,
        yminmax=yminmax,
        **kwargs
    )

    changes, pvals = zip(*[
        (x, y)
        for x, y in zip(changes, pvals)
        if x > xminmax[0] and x < xminmax[1] and
        y > yminmax[0] and y < yminmax[1]
    ])

    ax.scatter(changes, pvals, c="lightblue", zorder=0, alpha=0.3)

    f.savefig(
        os.path.join(
            data.name,
            re.sub(r"[ ></\?]", "_", kwargs.get("title", "Filtered")) +
            "_Volcano.png",
        ),
        bbox_inches="tight", dpi=300,
        transparent=True,
    )
    return f, ax
