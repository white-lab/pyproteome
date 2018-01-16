
from collections import OrderedDict
import os
from unittest import TestCase
import itertools

import pylab

from pyproteome import (
    analysis, cluster, correlation, data_sets, logo, motif, motifs, paths,
    tables, volcano, phosphosite,
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
            dirname=paths.MS_SEARCHED_DIR,
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
        self.data["merge"] = data_sets.merge_data(
            [
                self.data["CKH1"],
                self.data["CKX2"],
                self.data["CKC1"],
            ],
            name="CK-p25",
        )

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

    def test_filter(self):
        for data in self.data.values():
            for f in [
                {"p": .1},
                {"fold": 1.5},
                {"fold": 1/1.5},
                {"asym_fold": 1.5},
                {"asym_fold": 1/1.5},
                {"sequence": "AVDSLVPIGR"},
                {"sequence": ["AVDSLVPIGR"]},
                {"protein": "Pkm"},
                {"protein": ["Pkm"]},
                {"ion_score": 10},
                {"isolation": 15},
                {"fn": lambda x: len(x["Sequence"]) > 5},
                {"series": data["Isolation Interference"] < data["Ion Score"]},
                {"missed_cleavage": 0},
                {"median_quant": 10000},
                {"only_validated": False},
                {"mod_types": [(None, "Phospho")]},
                {"p": .1, "inverse": True},
                {"motif": motif.Motif(".......xP......")},
                {"confidence": "Low"},
                {"confidence": "Medium"},
                {"confidence": "High"},
                {"group_a": "CK-p25 Hip", "group_b": "CK Hip"},
            ]:
                data.filter(f)

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
        for f in [
            dict(sequence="AVDSLVPIGR"),
            dict(protein="Pkm"),
        ]:
            for cmp_groups in [self.cmp_groups, None]:
                for data in self.data.values():
                    analysis.plot.plot_all(
                        data,
                        f=f,
                        individual=True,
                        between=True,
                        cmp_groups=cmp_groups,
                    )

                    analysis.plot.plot_together(
                        data,
                        f=f,
                        cmp_groups=cmp_groups,
                    )

    def test_correlate_data(self):
        for d1, d2 in itertools.combinations(self.data.values(), 2):
            correlation.correlate_data_sets(
                d1, d2,
            )

    def test_correlate_signal(self):
        for data in self.data.values():
            correlation.correlate_signal(
                data,
                data.dropna(how="any").psms
                .iloc[0][list(data.channels.values())],
            )

    def test_clustermap(self):
        for data in self.data.values():
            cluster.plot.hierarchical_heatmap(
                data,
            )

    def test_find_tfs(self):
        for data in self.data.values():
            analysis.pathways.find_tfs(data)

    def test_changes_table(self):
        for data in self.data.values():
            tables.changes_table(data)

    def test_volcano_plot(self):
        for data in self.data.values():
            volcano.volcano_plot(data)

    def test_plot_volcano_filtered(self):
        for data in self.data.values():
            volcano.plot_volcano_filtered(data, {"asym_fold": 1.001})

    def test_write_full_tables(self):
        tables.write_full_tables(
            self.data.values(),
        )

    def test_logo(self):
        for data in self.data.values():
            logo.make_logo(data, {"asym_fold": 1.001})

    def test_neighborhood(self):
        for data in self.data.values():
            motifs.neighborhood.enriched_neighborhood(
                data,
                {"asym_fold": 1.001},
                "ACDEFGHIKLMNPQRSTVWY",
            )

    # def test_plogo(self):
    #     motifs.plogo.make_logo(self.data["CKH1"], {"asym_fold": 1.05})

    def test_icelogo(self):
        for data in self.data.values():
            motifs.icelogo.make_logo(data, {"asym_fold": 1.001})

    # def test_weblogo(self):
    #     for data in self.data.values():
    #         motifs.weblogo.make_logo(data)

    def test_phosphosite_enriched(self):
        for data in self.data.values():
            phosphosite.enriched(
                data,
                species="Mouse",
            )

    def test_cluster(self):
        for d in self.data.values():
            data = cluster.get_data(d)

            cluster.plot.pca(data)

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
        for data in self.data.values():
            cluster.auto.auto_clusterer(
                data,
                cluster_kwargs={"n_clusters": 10},
            )
