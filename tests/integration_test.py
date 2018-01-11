
from collections import OrderedDict
import os
import requests
from unittest import TestCase

from pyproteome import (
    analysis, data_sets, levels, logo, paths, tables, volcano,
)


DATA_URL = "https://github.com/white-lab/pyproteome-data/raw/master/test_data/"

DATAS = {
    "CKH1": "CKH1-pY-sup.msf",
    "CKX2": "CKX2-pY-sup.msf",
    "CKC1": "CKC1-pY-sup.msf",
}


class IntegrationTest(TestCase):
    def _fetch_data(self):
        try:
            os.makedirs(paths.MS_SEARCHED_DIR)
        except:
            pass

        for _, filename in DATAS.items():
            url = DATA_URL + filename
            response = requests.get(url, stream=True)
            response.raise_for_status()

            out_path = os.path.join(paths.MS_SEARCHED_DIR, filename)

            with open(out_path, mode="wb") as f:
                for block in response.iter_content(1024):
                    f.write(block)

    def setUp(self):
        paths.set_base_dir(os.path.abspath("."))
        self._fetch_data()

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

        self.data = {
            name: data_sets.DataSet(
                mascot_name=os.path.splitext(filename)[0],
                name=name,
                channels=self.channels[name],
                groups=self.groups,
            )
            for name, filename in DATAS.items()
        }

    def test_normalize_data(self):
        self.levels = {
            name: levels.get_channel_levels(data)
            for name, data in self.data.items()
        }

        self.data = {
            name: data.normalize(self.levels[name])
            for data, name in self.data.items()
        }

    def test_merge_data(self):
        self.test_normalize_data()

        merge = data_sets.merge_data(
            [
                self.data["CKH1"],
                self.data["CKX2"],
                self.data["CKC1"],
            ],
            name="CK-p25",
        )

        print(merge.shape[0])

        for data in [
            self.data["CKH1"],
            self.data["CKX2"],
            self.data["CKC1"],
            merge,
        ]:
            data.print_stats()

        # self.assertEqual(
        #     merge.shape[0],
        #     100,
        # )
        return merge

    def test_correlate_data(self):
        self.test_normalize_data()

        analysis.correlate_data_set(
            self.data["CKH1"],
            self.data["CKX2"],
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
            volcano.plot_volcano_filtered(data, {"asym_fold_cutoff": 1.25})

    def test_write_full_tables(self):
        merge = self.test_merge_data()

        tables.write_full_tables(
            [data for data in self.data.values()] + [merge]
        )

    def test_logo(self):
        self.test_normalize_data()

        for _, data in self.data.items():
            logo.make_logo(data, {"asym_fold_cutoff": 1.25})
