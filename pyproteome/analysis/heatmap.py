
import os

import numpy as np
import seaborn as sns

import pyproteome as pyp


def hierarchical_heatmap(
    data,
    cmp_groups=None,
    baseline_channels=None,
    minmax=None,
    show_y=False,
    title=None,
    filename="Hierarchical Heatmap.png",
    folder_name=None,
    **kwargs
):
    """
    Plot a hierarhically-clustered heatmap of a data set.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.DataSet`
    baseline_channels : list of str, optional
        List of channels to average and use as baseline for each row.
    row_cluster : bool, optional
    col_cluster : bool, optional

    Returns
    -------
    map : :class:`seaborn.ClusterGrid`
    """
    folder_name = pyp.utils.make_folder(
        data=data,
        folder_name=folder_name,
        sub="Clusters",
    )

    if cmp_groups is None:
        cmp_groups = [list(data.groups.keys())]
        # cmp_groups = [data.group_a, data.group_b]

    flat_cmp_groups = pyp.utils.flatten_list(cmp_groups)

    colors = [
        'red',
        'orange',
        'purple',
        'blue',
        'cyan',
        'green',
    ]

    group_colors = {
        data.channels[channel]: colors[
            flat_cmp_groups.index(group) % len(colors)
        ]
        for groups in cmp_groups
        for group in groups
        for channel in data.groups[group]
        if channel in data.channels
    }
    channels = list(group_colors.keys())

    raw = data.psms
    raw.index = data.psms.apply(
        lambda x:
        "{0}{1}{2}{3}{3}".format(
            " / ".join(x["Proteins"].genes),
            " : "
            if len(x['Modifications'].get_mods(['S', 'T', 'Y', 'M'])) > 0 else
            "",
            x["Modifications"].get_mods(['S', 'T', 'Y', 'M']).__str__(
                prot_index=0,
            ),
            # [
            #     r"$\textbf{{{0}}}${1}$\textbf{{{2}}}${3}$\textbf{{{4}}}$".format(
            #         i[0],
            #         i[1],
            #         i[2],
            #         i[3:6],
            #         i[6],
            #     )
            #     for i in pyp.motifs.generate_n_mers(
            #         x['Sequence'],
            #         letter_mod_types=["S", "T"],
            #         all_matches=False,
            #         n=11,
            #     )
            # ][0],
            "",
        ),
        axis=1,
    )
    raw = raw[~raw.index.duplicated(keep='first')]
    # raw = raw.sort_values("Fold Change", ascending=False)
    raw = raw[channels]

    def zscore(x):
        return (x - np.nanmean(x)) / np.nanstd(x)

    # raw = raw.apply(zscore, axis=1)
    raw = raw.apply(np.log2, axis=1)
    raw = raw.dropna(how="all")
    raw = raw.T.dropna(how="all").T
    raw = raw.fillna(0)

    if minmax is None:
        minmax = max([
            abs(raw.min().min()),
            abs(raw.max().max()),
        ])
    else:
        minmax = np.log2(minmax)

    raw = raw.fillna(value=-minmax)

    if raw.shape[0] < 1:
        return

    map = sns.clustermap(
        raw,
        col_colors=[group_colors[i] for i in raw.columns],
        vmin=-minmax,
        vmax=minmax,
        robust=True,
        **kwargs
    )
    map.ax_heatmap.tick_params(
        # direction='out',
        axis='x',
        width=1.5,
        length=15,
        bottom=True,
        # colors='k',
    )

    map.ax_heatmap.set_xticklabels(
        map.ax_heatmap.get_xticklabels(),
        rotation=45,
        horizontalalignment="right",
    )

    if title is not None:
        map.fig.suptitle(title)

    if not show_y:
        map.ax_heatmap.set_yticklabels([])
        map.ax_heatmap.set_ylabel("")

    if filename:
        map.savefig(
            os.path.join(folder_name, filename),
            bbox_inches="tight",
            dpi=pyp.DEFAULT_DPI,
            transparent=True,
        )

    return map
