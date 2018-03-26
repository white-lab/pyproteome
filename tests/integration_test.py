
from collections import OrderedDict
import os
from unittest import TestCase
import itertools
import shutil

import pylab

import pyproteome as pyp
from pyproteome import (
    analysis, cluster, correlation, data_sets, logo, motif, motifs, paths,
    tables, volcano, phosphosite,
)


from . import utils


DATA_URL = "https://github.com/white-lab/pyproteome-data/raw/master/test_data/"

DATAS = (
    "CK-H1-Global",
    "CK-X2-Global",
    "CK-C1-Global",
)


class IntegrationTest(TestCase):
    @classmethod
    def setUpClass(cls):
        paths.set_base_dir(
            os.path.abspath(
                os.path.join(".", "tests_output")
            )
        )

        utils.fetch_data(
            dirname=paths.MS_SEARCHED_DIR,
            datas=[i + ".msf" for i in DATAS],
            base_url=DATA_URL,
        )

        pylab.rcParams['figure.max_open_warning'] = 0
        pyp.DEFAULT_DPI = 100

    def setUp(self):
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
            "CK-H1-Global": ckh_channels,
            "CK-X2-Global": ckx_channels,
            "CK-C1-Global": ckc_channels,
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

        self.data = data_sets.load_all_data(
            chan_mapping=self.channels,
            groups=self.groups,
            check_raw=False,
        )
        self.data["merge"] = data_sets.merge_data(
            [
                self.data["CK-H1-Global"],
                self.data["CK-X2-Global"],
                self.data["CK-C1-Global"],
            ],
            name="CK-p25",
        )

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(paths.FIGURES_DIR, ignore_errors=True)

    def test_columns(self):
        for _, data in self.data.items():
            for col in data_sets.DATA_SET_COLS:
                self.assertIn(col, data.psms.columns)

            for _, col in data.channels.items():
                self.assertIn(col, data.psms.columns)

    def test_dropna(self):
        self.data = {
            name: data_sets.DataSet(
                name=name,
                channels=self.channels[name],
                groups=self.groups,
                dropna=True,
                check_raw=False,
            )
            for name in DATAS
        }

    def test_merge_subsets(self):
        for data in self.data.values():
            data.merge_subsequences()

    def test_add_data(self):
        data = self.data["CK-H1-Global"]
        data += self.data["CK-X2-Global"]
        data += self.data["CK-C1-Global"]

    def test_normalize_data(self):
        self.data = {
            name: data.normalize(data)
            for name, data in self.data.items()
        }

    def test_cmp_groups(self):
        self.data = {
            name: data.norm_cmp_groups(
                [
                    ["CK-p25 Hip", "CK Hip"],
                    ["CK-p25 Cortex", "CK Cortex"],
                    ["CK-p25 Cere", "CK Cere"],
                ]
            )
            for name, data in self.data.items()
        }

    def test_filter(self):
        for data in self.data.values():
            for f in [
                {"ambiguous": True},
                {"ambiguous": False},
                {"p": .1},
                {"q": .1},
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
                self.data["CK-H1-Global"],
                self.data["CK-X2-Global"],
                self.data["CK-C1-Global"],
            ],
            name="CK-p25",
        )

        merge = data_sets.merge_data(
            [
                self.data["CK-H1-Global"],
                self.data["CK-X2-Global"],
                self.data["CK-C1-Global"],
            ],
            name="CK-p25",
        )

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
        self.test_normalize_data()
        self.data["merge"] = self.data["merge"].norm_cmp_groups(
            self.cmp_groups,
        )
        for data in self.data.values():
            correlation.correlate_signal(
                data,
                data.dropna(how="any").psms
                .iloc[0][list(data.channels.values())],
            )

    def test_gsea(self):
        for data in self.data.values():
            analysis.pathways.gsea(
                data,
                metric="zscore",
                n_cpus=1,
                min_hits=15,
                pval=True,
                p_iter=10,
            )

    def test_psea(self):
        for data in self.data.values():
            analysis.pathways.gsea(
                data,
                metric="zscore",
                n_cpus=1,
                min_hits=15,
                p_sites=True,
                pval=True,
                p_iter=10,
            )

    def test_clustermap(self):
        for data in self.data.values():
            cluster.plot.hierarchical_heatmap(
                data,
            )

    def test_changes_table(self):
        for data in self.data.values():
            tables.changes_table(data)

    def test_plot_volcano(self):
        for data in self.data.values():
            volcano.plot_volcano(data)

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

    # def test_icelogo(self):
    #     for data in self.data.values():
    #         motifs.icelogo.make_logo(data, {"asym_fold": 1.001})

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

            cluster.plot.cluster_range(data, max_clusters=5)

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
                protein="Mbp",
            )

            cluster.plot.show_peptide_clusters(
                data, y_pred,
                [
                    {"seq": "AAVPDAVGK"},
                    {"protein": "Mbp"},
                ],
            )

    def test_auto_cluster(self):
        for data in self.data.values():
            for d in [
                data,
                data.norm_cmp_groups(self.cmp_groups),
            ]:
                cluster.auto.auto_clusterer(
                    d,
                    cluster_kwargs={"n_clusters": 5},
                )
