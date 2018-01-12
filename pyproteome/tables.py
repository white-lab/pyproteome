
from __future__ import division

import os
import re

import pandas as pd

from . import utils
from .motifs import motif


def _prep_csv(data, postfix="table", folder_name=None, csv_name=None):
    data_name = getattr(data, "name", "Data")
    if folder_name is None:
        folder_name = data_name

    if csv_name is None:
        csv_name = "{}-{}-{}.csv".format(
            re.sub("[ ></]", "_", data_name),
            re.sub("[ ></]", "_", data.enrichment),
            postfix,
        )

    utils.make_folder(folder_name)

    return os.path.join(folder_name, csv_name)


def _get_table_title(f=None, running_title=None):
    if running_title is None:
        running_title = []

    if f is not None:
        if 'asym_fold_cutoff' in f:
            running_title.append(
                "Upregulated" if f["asym_fold_cutoff"] > 1 else "Downregulated"
            )

        if "p_cutoff" in f:
            running_title.append(
                "p-{:.3e}".format(f["p_cutoff"])
            )

        if "group_a" in f or "group_b" in f:
            running_title.append(
                "{}vs{}".format(
                    f.get("group_a", ""),
                    f.get("group_b", ""),
                )
            )

    return "-".join(running_title)


def motif_table(
    data, f,
    p_cutoff=0.05,
    sort="p-value",
    folder_name=None, csv_name=None,
    **kwargs
):
    csv_name = _prep_csv(
        data,
        folder_name=folder_name,
        csv_name=csv_name,
        postfix=_get_table_title(f=f, running_title=["motifs"]),
    )

    hits = motif.run_motif_enrichment(
        data, f,
        **kwargs
    )[0]

    hits = hits.sort_values(
        sort,
        ascending=True if sort == "p-value" else False,
    )

    if csv_name:
        hits.to_csv(csv_name)

    hits = hits[
        hits[
            "pp-value" if kwargs.get("pp_value", False) else "p-value"
        ] < p_cutoff
    ]

    return hits.style.set_table_styles([
        {"selector": "*", "props": [("font-family", "monospace")]},
        {"selector": "th:first-child", "props": [("display", "none")]},
    ])


def changes_table(
    data,
    sort="p-value",
    folder_name=None, csv_name=None,
):
    """
    Show a table of fold changes and p-values.

    Parameters
    ----------
    data : :class:`DataSet<pyproteome.data_sets.DataSet>`
    sort : str, optional
    folder_name : str, optional
    csv_name : str, optional
    """
    csv_name = _prep_csv(
        data,
        folder_name=folder_name,
        csv_name=csv_name,
        postfix=_get_table_title(running_title=["changes"]),
    )
    psms = getattr(data, "psms", data)

    psms = psms[
        [
            "Proteins", "Sequence", "Modifications",
            "Fold Change", "p-value", "Validated",
        ]
    ].copy()
    psms["Sequence"] = psms["Sequence"].apply(
        lambda x: "{} ({})".format(x, x.modifications)
    )

    if sort == "Fold Change":
        psms["Fold Change-Sort"] = psms["Fold Change"].apply(
            lambda x: max([x, 1 / x])
        )
        psms.sort_values("Fold Change-Sort", inplace=True, ascending=False)
        psms.drop("Fold Change-Sort", axis=1, inplace=True)
    else:
        psms.sort_values(sort, inplace=True, ascending=True)

    psms.drop("Modifications", axis=1, inplace=True)

    if csv_name:
        psms.to_csv(csv_name)

    # back_colors = {
    #     True: "#BBFFBB",  # light green
    #     False: "#FFBBBB",  # light red
    # }

    if psms.empty:
        return psms

    # return psms.style.apply(  # Color validated rows
    #     lambda row: [
    #         "background-color: " + back_colors[row["Validated"]]
    #         for _ in row
    #     ],
    #     axis=1,
    # )
    return psms.style.set_table_styles(  # Hide index and "Validated" columns
        [
            {"selector": "th:first-child", "props": [("display", "none")]},
            {"selector": "td:last-child", "props": [("display", "none")]},
            {"selector": "th:last-child", "props": [("display", "none")]},
            {"selector": "*", "props": [("text-align", "left")]},
        ]
    )


def _rewrite_enrichments(lst):
    mapping = {
        "MPM2": "pST",
        "MAPK / CDK": "pST",
        "MAPK + CDK": "pST",
        "NTA": "pST",
        "IMAC": "pST",
    }
    return "+".join(
        sorted(
            set(
                mapping.get(i, i)
                for i in lst
            )
        )
    )


def write_full_tables(datas, folder_name="All", out_name="Full Data.xlsx"):
    """
    Write a full list of data sets to a single .xlsx file.

    Parameters
    ----------
    datas : list of :class:`DataSet<pyproteome.data_sets.DataSet>`
    folder_name : str, optional
    out_name : str, optional
    """

    utils.make_folder(folder_name)

    if folder_name is not None:
        out_name = os.path.join(folder_name, out_name)

    writer = pd.ExcelWriter(out_name, engine="xlsxwriter")

    for data in datas:
        channels = [
            (name, data.channels[name])
            for group in data.groups.values()
            for name in group
            if name in data.channels and
            data.channels[name] in data.psms.columns
        ]
        df = data.psms[
            [
                "Proteins", "Sequence",
                "Fold Change", "p-value", "Validated",
            ] + [
                chan
                for _, chan in channels
            ]
        ].copy()

        df.rename(
            columns={
                chan: name
                for name, chan in channels
            },
            inplace=True,
        )
        df.insert(
            2, "Modifications",
            df["Sequence"].apply(
                lambda x: str(x.modifications)
            ),
        )
        df["Sequence"] = df["Sequence"].apply(str)
        df.sort_values("p-value", inplace=True, ascending=True)

        ws_name = re.sub(
            "/",
            "+",
            "{} - {}".format(
                data.name,
                _rewrite_enrichments(data.enrichments),
            ),
        )
        df.to_excel(
            writer,
            sheet_name=ws_name,
            index=False,
        )

        ws = writer.sheets[ws_name]
        ws.freeze_panes(1, 0)
        ws.set_column(0, 0, 60)
        ws.set_column(1, 1, 30)
        ws.set_column(2, 2, 20)
        ws.set_column(3, 3, 12)
        ws.set_column(4, 4, 12)

    writer.save()
