
from collections import OrderedDict
import os
from unittest import TestCase

import pylab

from pyproteome import (
    analysis, cluster, correlation, data_sets, logo, paths, plot,
    tables, volcano, plogo, icelogo, weblogo, phosphosite,
)

from . import utils


DATA_URL = "https://github.com/white-lab/pyproteome-data/raw/master/test_data/"

DATAS = {
    "CKH1": "CKH1-pY-sup.msf",
    "CKX2": "CKX2-pY-sup.msf",
    "CKC1": "CKC1-pY-sup.msf",
}


class IntegrationTest(TestCase):

    def setUp(self):
        paths.set_base_dir(os.path.abspath("."))

        utils.fetch_data(
            dir=paths.MS_SEARCHED_DIR,
            datas=DATAS,
            base_url=DATA_URL,
        )

        pylab.rcParams['figure.max_open_warning'] = 0

        ck_channels = OrderedDict(
            [
                ("3130 CK",     "126"),
                ("3131 CK-p25", "127"),
                ("3145 CK-p25", "128"),
                ("3146 CK-p25", "129"),
                ("3148 CK",     "130"),
                ("3157 CK",     "131"),
            ]
        )

        ckh_channels = OrderedDict([
            ("{} Hip".format(key), val)
            for key, val in ck_channels.items()
        ])
        ckx_channels = OrderedDict([
            ("{} Cortex".format(key), val)
            for key, val in ck_channels.items()
        ])
        ckc_channels = OrderedDict([
            ("{} Cere".format(key), val)
            for key, val in ck_channels.items()
        ])
        self.channels = {
            "CKH1": ckh_channels,
            "CKX2": ckx_channels,
            "CKC1": ckc_channels,
        }
        self.groups = OrderedDict(
            [
                (
                    "CK-p25",
                    [
                        "3131 CK-p25",
                        "3145 CK-p25",
                        "3146 CK-p25",
                    ],
                ),
                (
                    "CK",
                    [
                        "3130 CK",
                        "3148 CK",
                        "3157 CK",
                    ],
                ),
            ]
        )
        self.groups.update([
            (
                "{} {}".format(key, tissue),
                ["{} {}".format(i, tissue) for i in val],
            )
            for tissue in ["Hip", "Cortex", "Cere"]
            for key, val in self.groups.items()
        ])
        self.cmp_groups = [
            ["CK-p25 Hip", "CK Hip"],
            ["CK-p25 Cortex", "CK Cortex"],
            ["CK-p25 Cere", "CK Cere"],
        ]

        self.data = {
            name: data_sets.DataSet(
                mascot_name=os.path.splitext(filename)[0],
                name=name,
                channels=self.channels[name],
                groups=self.groups,
            )
            for name, filename in DATAS.items()
        }

    def test_dropna(self):
        self.data = {
            name: data_sets.DataSet(
                mascot_name=os.path.splitext(filename)[0],
                name=name,
                channels=self.channels[name],
                groups=self.groups,
                dropna=True,
            )
            for name, filename in DATAS.items()
        }

    def test_merge_subsets(self):
        self.data = {
            name: data_sets.DataSet(
                mascot_name=os.path.splitext(filename)[0],
                name=name,
                channels=self.channels[name],
                groups=self.groups,
                merge_subsets=True,
            )
            for name, filename in DATAS.items()
        }

    def test_add_data(self):
        data = self.data["CKH1"]
        data += self.data["CKX2"]
        data += self.data["CKC1"]

    def test_normalize_data(self):
        self.data = {
            name: data.normalize(data)
            for name, data in self.data.items()
        }

    def test_merge_data(self):
        self.test_normalize_data()

        merge = data_sets.merge_data([], name="Empty")

        merge = data_sets.merge_data(
            [
                self.data["CKH1"],
                self.data["CKX2"],
                self.data["CKC1"],
            ],
            name="CK-p25",
            merge_subsets=True,
        )

        merge = data_sets.merge_data(
            [
                self.data["CKH1"],
                self.data["CKX2"],
                self.data["CKC1"],
            ],
            name="CK-p25",
        )

        with open(os.devnull, 'w') as f:
            for data in [
                self.data["CKH1"],
                self.data["CKX2"],
                self.data["CKC1"],
                merge,
            ]:
                data.print_stats(out=f)

        # self.assertEqual(
        #     merge.shape[0],
        #     100,
        # )
        return merge

    def test_plot_all(self):
        merge = self.test_merge_data()

        for f in [
            dict(sequence="AVDSLVPIGR"),
            dict(protein="Pkm"),
        ]:
            for cmp_groups in [self.cmp_groups, None]:
                plot.plot_all(
                    merge,
                    f=f,
                    individual=True,
                    between=True,
                    cmp_groups=cmp_groups,
                )

                plot.plot_together(
                    merge,
                    f=f,
                    cmp_groups=cmp_groups,
                )

    def test_correlate_data(self):
        self.test_normalize_data()

        correlation.correlate_data_sets(
            self.data["CKH1"],
            self.data["CKX2"],
        )

        correlation.correlate_signal(
            self.data["CKH1"],
            self.data["CKH1"][self.data["CKH1"]["Sequence"] == "AVDSLVPIGR"],
        )

    def test_clustermap(self):
        self.test_normalize_data()

        for data in self.data.values():
            analysis.hierarchical_heatmap(
                data,
            )

    def test_changes_table(self):
        self.test_normalize_data()

        for _, data in self.data.items():
            tables.changes_table(data)

    def test_volcano_plot(self):
        self.test_normalize_data()

        for _, data in self.data.items():
            volcano.volcano_plot(data)

    def test_plot_volcano_filtered(self):
        self.test_normalize_data()

        for _, data in self.data.items():
            volcano.plot_volcano_filtered(data, {"asym_fold": 1.25})

    def test_write_full_tables(self):
        merge = self.test_merge_data()

        tables.write_full_tables(
            [data for data in self.data.values()] + [merge]
        )

    def test_logo(self):
        self.test_normalize_data()

        for _, data in self.data.items():
            logo.make_logo(data, {"asym_fold": 1.05})

    # def test_plogo(self):
    #     self.test_normalize_data()
    #
    #     plogo.make_logo(self.data["CKH1"], {"asym_fold": 1.05})

    def test_icelogo(self):
        self.test_normalize_data()

        icelogo.make_logo(self.data["CKH1"], {"asym_fold": 1.05})

    # def test_weblogo(self):
    #     self.test_normalize_data()
    #
    #     for _, data in self.data.items():
    #         weblogo.make_logo(data, {})

    def test_phosphosite_enriched(self):
        phosphosite.enriched(
            self.data["CKH1"],
            species="Mouse",
        )

    def test_cluster(self):
        merge = self.test_merge_data()

        data = cluster.get_data(merge)

        cluster.pca(data)

        _, y_pred = cluster.cluster(
            data,
            n_clusters=6,
        )

        cluster.cluster_range(data, max_clusters=5)

        cluster.plot.cluster_corrmap(data, y_pred)

        y_pred = cluster.cluster_clusters(data, y_pred)

        cluster.plot.cluster_corrmap(data, y_pred)

        f, ax = pylab.subplots()
        cluster.plot.cluster_corrmap(data, y_pred, ax=ax)

        cluster.plot.plot_all_clusters(data, y_pred)

        cluster.plot.plot_cluster(
            data, y_pred,
            cluster_n=0,
        )

        cluster.plot.show_cluster(
            data, y_pred,
            protein="Pkm",
        )

        cluster.plot.show_peptide_clusters(
            data, y_pred,
            [
                {"seq": "AVDSLVPIGR"},
                {"protein": "Pkm"},
            ],
        )

    def test_auto_cluster(self):
        merge = self.test_merge_data()

        cluster.auto.auto_clusterer(
            merge,
            cluster_kwargs={"n_clusters": 10},
        )
