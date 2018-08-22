# -*- coding: UTF-8 -*-
"""
This module provides functionality for data set analysis.

Functions include volcano plots, sorted tables, and plotting sequence levels.
"""

from __future__ import absolute_import, division

import logging
import os

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import pyproteome as pyp


LOGGER = logging.getLogger("pyproteome.correlation")


def correlate_data_sets(
    data1, data2,
    older_name=None,
    folder_name=None,
    filename=None,
    adjust=True,
    label_cutoff=1.5,
    show_labels=True,
):
    """
    Plot the correlation between peptides levels in two different data sets.

    Parameters
    ----------
    data1 : :class:`pyproteome.data_sets.DataSet`
    data2 : :class:`pyproteome.data_sets.DataSet`
    folder_name : str, optional
    filename : str, optional
    """
    folder_name = pyp.utils.make_folder(
        data=data1,
        folder_name=folder_name,
        sub="Correlation Analysis",
    )

    merged = pd.merge(
        data1.psms, data2.psms,
        on="Sequence",
    ).dropna(subset=("Fold Change_x", "Fold Change_y"))

    f, ax = plt.subplots(
        figsize=(4, 3),
    )

    sns.regplot(
        x=np.log2(merged["Fold Change_x"]),
        y=np.log2(merged["Fold Change_y"]),
        ax=ax,
        scatter_kws={
            "alpha": .5,
        }
    )

    label_cutoff = np.log2(label_cutoff)

    if show_labels:
        texts, x_s, y_s = [], [], []

        for index, row in merged.iterrows():
            x = row["Fold Change_x"]
            y = row["Fold Change_y"]
            ratio = np.log2(x / y)
            x, y = np.log2(x), np.log2(y)

            if ratio < label_cutoff and ratio > - label_cutoff:
                continue

            x_s.append(x)
            y_s.append(y)

            txt = " / ".join(row["Proteins_x"].genes)
            txt = txt[:20] + ("..." if len(txt) > 20 else "")

            text = ax.text(
                x, y, txt,
            )

            text.set_bbox(
                dict(
                    color="lightgreen" if ratio < 0 else "pink",
                    alpha=0.8,
                )
            )
            texts.append(text)

        if adjust:
            pyp.utils.adjust_text(
                x=x_s,
                y=y_s,
                texts=texts,
                ax=ax,
                lim=400,
                force_text=0.1,
                force_points=0.1,
                arrowprops=dict(arrowstyle="->", relpos=(0, 0), lw=1),
                only_move={
                    "points": "y",
                    "text": "xy",
                }
            )

    min_x = min(np.log2(merged["Fold Change_x"]))
    max_x = max(np.log2(merged["Fold Change_x"]))
    ax.plot([min_x, max_x], [min_x, max_x], "--", color="k")
    # ax.plot(
    #     [min_x, max_x],
    #     [min_x + label_cutoff, max_x + label_cutoff],
    #     color="lightgreen",
    #     linestyle=":",
    # )
    # ax.plot(
    #     [min_x, max_x],
    #     [min_x - label_cutoff, max_x - label_cutoff],
    #     color="pink",
    #     linestyle=":",
    # )

    name1 = data1.name
    name2 = data2.name

    ax.set_xlabel("Fold Change -- {}".format(name1))
    ax.set_ylabel("Fold Change -- {}".format(name2))

    ax.set_xticklabels(
        ["{:.2f}".format(i) for i in np.power(2, ax.get_xticks())],
        # fontsize=20,
    )
    ax.set_yticklabels(
        ["{:.2f}".format(i) for i in np.power(2, ax.get_yticks())],
        # fontsize=20,
    )

    pear_corr = merged["Fold Change_x"].corr(
        merged["Fold Change_y"],
        method="pearson",
    )
    # spear_corr = merged["Fold Change_x"].corr(
    #     merged["Fold Change_y"],
    #     method="spearman",
    # )

    ax.set_title(
        (
            r"Pearson's: $\rho$={:.2f}"
            # r"Spearman's: $\rho$={:.2f}"
        ).format(pear_corr)
    )

    if filename:
        f.savefig(
            os.path.join(folder_name, filename),
            transparent=True,
            dpi=pyp.DEFAULT_DPI,
            bbox_inches="tight",
        )


def _scatter_plots(
    cp, signal,
    data_chans, signal_chans, signal_groups,
    scatter_cols=4,
    xlabel="",
    scatter_colors=None,
    scatter_symbols=None,
):
    scatter_colors = scatter_colors or {}
    scatter_symbols = scatter_symbols or {}

    if cp.shape[0] < 1:
        return None

    df = pd.DataFrame(
        [
            (
                row,
                str(row["Sequence"]),
                signal[sig_chan],
                row[data_chan],
                sig_group,
                scatter_colors.get(sig_group, "black"),
            )
            for _, row in cp.psms.sort_values("Correlation").iterrows()
            for data_chan, sig_chan, sig_group in zip(
                data_chans, signal_chans, signal_groups,
            )
        ],
        columns=("row", "sequence", "x", "y", "group", "color"),
    )
    g = sns.lmplot(
        x="x",
        y="y",
        hue="group",
        col="sequence",
        col_wrap=scatter_cols,
        data=df,
        sharey=False,
        fit_reg=False,
    )
    for ax in g.axes:
        seq = ax.get_title().split("=", 1)[1].strip()
        df_cp = df[df["sequence"] == seq]
        row = df_cp["row"].iloc[0]

        row_title = " / ".join(str(i.gene) for i in row["Proteins"])
        row_seq = str(row["Sequence"])
        row_mods = str(row["Modifications"].get_mods([(None, "Phospho")]))

        sns.regplot(
            x="x",
            y="y",
            data=df_cp,
            scatter=False,
            ax=ax,
        )

        ax.set_title(
            row_title[:20] + ("..." if len(row_title) > 20 else ""),
            # fontsize=28,
            # fontweight="bold",
        )
        # ax.set_xlabel(
        #     "{}\n$\\rho = {:.2f}$".format(
        #         xlabel,
        #         row["Correlation"],
        #     ),
        #     fontsize=22,
        # )
        ax.set_ylabel(
            row_seq[:20] + (
                "..." if len(row_seq) > 20 else ""
            ) + (
                "\n({})".format(
                    row_mods[:20] + ("..." if len(row_mods) > 20 else "")
                ) if row_mods else ""
            ),
            # fontsize=20,
        )

        # for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
        #     tick.label.set_fontsize(20)

        ax.legend(loc="best")

    return g


def _remove_lesser_dups(labels, compress_sym=False):
    new_labels = []
    labels = list(labels)

    for index, (x, y, label) in enumerate(labels):
        for o_index, (o_x, o_y, o_label) in enumerate(labels):
            if index == o_index:
                continue
            if label != o_label:
                continue
            if not compress_sym and (y < 0) != (y < 0):
                continue

            mul = 1 if y >= 0 else -1
            o_mul = 1 if o_y >= 0 else -1

            if (
                mul * y < o_mul * o_y
            ) or (
                (
                    mul * y <= o_mul * o_y
                ) and
                index < o_index
            ):
                break
        else:
            new_labels.append((x, y, label))

    return new_labels


def correlate_signal(
    data, signal,
    corr_cutoff=0.8,
    scatter_cols=4,
    options=None,
    folder_name=None,
    title=None,
    show_duplicates=False,
    scatter_colors=None,
    scatter_symbols=None,
    show_scatter=True,
    figsize=(12, 10),
    xlabel="",
):
    """
    Calculate the correlation between levels of each peptide in a data set and
    a given signal variable.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.DataSet`
    signal : :class:`pandas.Series`

    Returns
    -------
    f_corr : :class:`matplotlib.figure.Figure`
    f_scatter : :class:`matplotlib.figure.Figure`
    """
    folder_name = pyp.utils.make_folder(
        data=data,
        folder_name=folder_name,
        sub="Correlation Analysis",
    )

    options = options or {}

    highlight = options.get('highlight', {})
    hide = options.get('hide', {})
    edgecolors = options.get('edgecolors', {})
    rename = options.get('rename', {})
    scatter_colors = options.get('scatter_colors', {})
    scatter_symbols = options.get('scatter_symbols', {})

    if title is None:
        title = data.name

    cp = data.copy()
    cp = cp.rename_channels()

    cols = getattr(signal, "columns", getattr(signal, "index", []))

    signal_groups = [
        label
        for label, group in cp.groups.items()
        for chan in group
        if chan in cp.channels.keys() and
        chan in cols
    ]

    signal_chans = [
        chan
        for group in cp.groups.values()
        for chan in group
        if chan in cp.channels.keys() and
        chan in cols
    ]
    data_chans = [
        cp.channels[chan]
        for chan in signal_chans
    ]

    cp.psms["Correlation"] = cp.psms.apply(
        lambda row:
        signal[signal_chans].corr(
            pd.to_numeric(row[data_chans]),
            method="spearman",
            min_periods=5,
        ),
        axis=1,
    )

    f_corr, ax = plt.subplots(figsize=figsize)
    x, y, colors = [], [], []
    sig_x, sig_y, sig_labels = [], [], []

    for index, (_, row) in enumerate(cp.psms.iterrows()):
        if (
            row["Correlation"] == 0 or
            np.isinf(row["Correlation"]) or
            np.isnan(row["Correlation"])
        ):
            continue

        x.append(index)
        y.append(row["Correlation"])

        sig = abs(row["Correlation"]) >= corr_cutoff

        if sig:
            sig_x.append(index)
            sig_y.append(row["Correlation"])
            sig_labels.append(" / ".join(sorted(row["Proteins"].genes)))

        colors.append(
            "blue"
            if sig else
            "{:.2f}".format(
                max([len(row[data_chans].dropna()) / len(data_chans) - .25, 0])
            )
        )

    ax.scatter(x, y, c=colors)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    labels = zip(sig_x, sig_y, sig_labels)
    LOGGER.info("Showing names for {} genes".format(len(sig_x)))

    if not show_duplicates:
        labels = _remove_lesser_dups(labels)

    texts = []

    for xs, ys, txt in labels:
        if any(i.strip() in hide for i in txt.split("/")):
            continue

        txt = rename.get(txt, txt)

        text = ax.text(
            xs, ys,
            txt[:20] + ("..." if len(txt) > 20 else ""),
        )

        if txt in highlight:
            text.set_fontsize(20)

        text.set_bbox(
            dict(
                facecolor="lightgreen" if xs > 0 else "pink",
                alpha=1,
                linewidth=0.5 if txt not in edgecolors else 3,
                edgecolor=edgecolors.get(txt, "black"),
                boxstyle="round",
            )
        )

        texts.append(text)

    pyp.utils.adjust_text(
        x=sig_x,
        y=sig_y,
        texts=texts,
        ax=ax,
        lim=400,
        force_text=0.3,
        force_points=0.01,
        arrowprops=dict(arrowstyle="->", relpos=(0, 0), lw=1),
        only_move={
            "points": "y",
            "text": "xy",
        }
    )

    ax.set_xlabel(
        "Index",
        fontsize=20,
    )
    # ax.set_yticklabels(
    #     "{:.3}".format(i)
    #     for i in np.power(1/10, ax.get_yticks())
    # )
    ax.set_ylabel(
        "Correlation",
        fontsize=20,
    )

    if title:
        ax.set_title(title, fontsize=32)

    cp.psms = cp.psms[cp.psms["Correlation"].apply(abs) >= corr_cutoff]

    if f_corr:
        f_corr.savefig(
            os.path.join(folder_name, "Correlation.png"),
            bbox_inches="tight",
            dpi=pyp.DEFAULT_DPI,
            transparent=True,
        )

    f_scatter = None

    if show_scatter:
        f_scatter = _scatter_plots(
            cp, signal, data_chans, signal_chans, signal_groups,
            scatter_cols=scatter_cols,
            scatter_colors=scatter_colors,
            scatter_symbols=scatter_symbols,
        )

        if f_scatter:
            f_scatter.savefig(
                os.path.join(folder_name, "Correlation Scatter.png"),
                bbox_inches="tight",
                # dpi=pyp.DEFAULT_DPI,
                dpi=pyp.DEFAULT_DPI / 6,
                transparent=True,
            )

    return f_corr, f_scatter
