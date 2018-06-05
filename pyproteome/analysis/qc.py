
from collections import OrderedDict
import itertools

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


def _is_qc_equal(row_a, row_b):
    if row_a["Proteins"].genes != row_b["Proteins"].genes:
        return False

    str_a = str(row_a["Sequence"])
    str_b = str(row_b["Sequence"])

    if str_a == str_b:
        return False

    if str_a not in str_b and str_b not in str_a:
        return False

    try:
        ind_a = str_a.index(str_b)
        slice_a = str_a[:ind_a] + str_a[ind_a + len(str_b):]
    except ValueError:
        slice_a = ""

    try:
        ind_b = str_b.index(str_a)
        slice_b = str_b[:ind_b] + str_b[ind_b + len(str_a):]
    except ValueError:
        slice_b = ""

    if any([i.islower() for i in slice_a + slice_b]):
        return False

    return True


def peptide_corr(ds):
    ds = ds.filter(fn=lambda x: len(x["Proteins"]) < 2)

    strs = set(ds["Sequence"].apply(str))

    def _fltr(x):
        x_str = str(x["Sequence"])
        return bool(any([i in x_str for i in strs if i != x_str]))

    print(ds.shape)
    ds = ds.filter(fn=_fltr)
    print(ds.shape)

    corrs = []

    chans = [
        val
        for key, val in ds.channels.items()
        if any([key in grp for grp in ds.groups.values()])
    ]
    weight_chans = [
        key + "_weight" if (key + "_weight") in ds.psms.columns else val
        for key, val in ds.channels.items()
        if any([key in grp for grp in ds.groups.values()])
    ]

    for (ind_a, row_a), (ind_b, row_b) in itertools.combinations(
        ds.psms.iterrows(), 2,
    ):
        if ind_a == ind_b:
            continue

        if not _is_qc_equal(row_a, row_b):
            continue

        data_a = pd.to_numeric(row_a[chans])
        data_b = pd.to_numeric(row_b[chans])

        weights_a = pd.to_numeric(row_a[weight_chans])
        weights_b = pd.to_numeric(row_b[weight_chans])

        s_corr = data_a.corr(
            data_b,
            method="spearman",
            min_periods=5,
        )
        p_corr = data_a.corr(
            data_b,
            method="pearson",
            min_periods=5,
        )

        if np.isnan(s_corr) or np.isnan(p_corr):
            continue

        corrs.append(
            OrderedDict([
                ("protein", row_a["Proteins"]),
                ("seq_a", row_a["Sequence"]),
                ("seq_b", row_b["Sequence"]),
                ("spearman", s_corr),
                ("pearson", p_corr),
                (
                    "score",
                    np.max([
                        row_a["Ion Score"],
                        row_b["Ion Score"],
                    ]),
                ),
                (
                    "isolation",
                    np.max([
                        row_a["Isolation Interference"],
                        row_b["Isolation Interference"],
                    ]),
                ),
                (
                    "quant",
                    np.log10(
                        np.mean([
                            np.nanmedian(weights_a),
                            np.nanmedian(weights_b),
                        ])
                    ),
                ),
                (
                    "mmc",
                    np.random.random() / 5 + np.max([
                        row_a["Missed Cleavages"],
                        row_b["Missed Cleavages"],
                    ]),
                ),
                (
                    "len",
                    np.max([
                        len(row_a["Sequence"].pep_seq),
                        len(row_b["Sequence"].pep_seq),
                    ]),
                ),
                (
                    "under",
                    np.random.random() / 5 + int(
                        row_a["Sequence"].is_underlabeled or
                        row_b["Sequence"].is_underlabeled
                    ),
                ),
            ])
        )

    df = pd.DataFrame(corrs).sort_values("pearson")

    if df.shape[0] > 10:
        f, ax = plt.subplots()
        sns.kdeplot(
            [i["spearman"] for i in corrs],
            shade=True,
        )

        sns.kdeplot(
            [i["pearson"] for i in corrs],
            shade=True,
        )
        ax.set_title(ds.name)

        for y in [
            "quant",
            "score",
            "isolation",
            "mmc",
            "sites",
            "len",
            "under",
        ]:
            sns.jointplot(
                x="pearson",
                y=y,
                data=df,
                alpha=10 / max([10, df.shape[0]]),
            ).plot_joint(sns.kdeplot, zorder=0, shade=True)

    return df
