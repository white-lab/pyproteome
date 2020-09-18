
from collections import OrderedDict

import numpy as np
import seaborn as sns
import re

import pyproteome as pyp


def _zscore(x):
    return (x - np.nanmean(x)) / np.nanstd(x)


def hierarchical_heatmap(
    data,
    cmp_groups=None,
    minmax=0,
    zscore=False,
    show_y=False,
    title=None,
    **kwargs
):
    """
    Plot a hierarhically-clustered heatmap of a data set.

    Parameters
    ----------
    data : :class:`pyproteome.data_sets.DataSet`
    cmp_groups : list of list of str
    minmax : float, optional
    zscore : bool, optional
    show_y : bool, optional
    title : str, optional
    kwargs : dict
        Kwargs passed directly to :func:`seaborn.clustermap`.

    Returns
    -------
    map : :class:`seaborn.ClusterGrid`
    """
    data = data.copy()

    if cmp_groups is None:
        # cmp_groups = [list(data.groups.keys())]
        cmp_groups = (
            data.cmp_groups or
            [i for i in [data.group_b, data.group_a] if i] or
            [list(data.groups.keys())]
        )

    flat_cmp_groups = pyp.utils.flatten_list(cmp_groups)

    colors = sns.color_palette(
        "hls", len(flat_cmp_groups),
    ).as_hex()

    group_colors = OrderedDict([
        (
            data.channels[channel],
            colors[
                flat_cmp_groups.index(group)
            ],
        )
        for groups in cmp_groups
        for group in [i for i in groups if i in data.groups]
        for channel in data.groups[group]
        if channel in data.channels
    ])
    channels = list(group_colors.keys())

    raw = data.psms
    raw.index = data.psms.apply(
        lambda x:
        "{0}{1}{2}{3}".format(
            pyp.utils.get_name(x["Proteins"]),
            " : "
            if len(x['Modifications'].get_mods(['S', 'T', 'Y', 'M'])) > 0 else
            "",
            re.sub(
                r'(\d+)',
                # r'$_{\1}$',
                r'\1',
                x["Modifications"].get_mods(['S', 'T', 'Y', 'M']).__str__(
                    prot_index=0,
                    show_mod_type=False,
                ),
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
            # ': ' + x["Sequence"].pep_seq,
        ),
        axis=1,
    )
    # raw = raw.sort_values("Fold Change", ascending=False)
    raw = raw[channels]
    sort_index = list(raw.index.drop_duplicates())
    raw = raw.groupby(raw.index).agg('median')
    raw['Sort Index'] = raw.index.map(lambda x: sort_index.index(x)) 
    raw = raw.sort_values('Sort Index')
    del raw['Sort Index']

    raw = raw.apply(np.log2, axis=1)
    if zscore:
        raw = raw.apply(_zscore, axis=1)

    raw = raw.replace([np.inf, -np.inf], np.nan)
    raw = raw.dropna(how="all")
    raw = raw.T.dropna(how="all").T

    if minmax == 0:
        minmax = max([
            abs(raw.min().min()),
            abs(raw.max().max()),
        ])
    else:
        minmax = np.log2(minmax)

    cluster_rowcol = (
        kwargs.get('row_cluster', True) or
        kwargs.get('col_cluster', True)
    )

    if cluster_rowcol:
        raw_na = raw.fillna(value=0)
    else:
        raw_na = raw.copy()

    if raw.shape[0] < 1:
        return

    row_colors = kwargs.pop('row_colors', None)
    col_colors = kwargs.pop('col_colors', None)

    if row_colors:
        row_colors = raw_na.index.map(row_colors).tolist()
        row_colors = np.array(row_colors).T.tolist()
    if col_colors is None:
        col_colors = (
            [group_colors[i] for i in raw.columns]
            if len(flat_cmp_groups) > 0 else None
        )
    elif callable(col_colors):
        col_colors = col_colors(raw.columns)

    if 'cmap' not in kwargs:
        kwargs['cmap'] = 'coolwarm'

    map = sns.clustermap(
        raw_na,
        row_colors=row_colors,
        col_colors=col_colors,
        vmin=-minmax,
        vmax=minmax,
        robust=True,
        **kwargs
    )

    kwargs.pop('metric', None)
    kwargs.pop('method', None)
    kwargs.pop('row_cluster', None)
    kwargs.pop('col_cluster', None)
    kwargs.pop('colors_ratio', None)
    kwargs.pop('figsize', None)
    kwargs.pop('cmap', None)

    if cluster_rowcol:
        raw_blank = raw.copy()

        if map.dendrogram_row:
            raw_blank = raw_blank.iloc[map.dendrogram_row.reordered_ind]

        if map.dendrogram_col:
            raw_blank = raw_blank[[
                raw_blank.columns[i]
                for i in map.dendrogram_col.reordered_ind
            ]]

        mask = raw_blank.isnull()
        raw_blank = raw_blank.fillna(0)

        sns.heatmap(
            raw_blank,
            mask=~mask,
            vmin=0,
            vmax=1,
            ax=map.ax_heatmap,
            cbar=False,
            cmap='binary',
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

    return map
