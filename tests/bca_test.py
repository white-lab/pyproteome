
import os
import requests
from unittest import TestCase

import pylab

from pyproteome import bca, paths


DATA_URL = "https://github.com/white-lab/pyproteome-data/raw/master/test_data/"

DATAS = {
    "CK-p25": "CK-p25-bca.xlsx",
}


class BCATest(TestCase):
    def _fetch_data(self):
        try:
            os.makedirs(paths.BCA_ASSAY_DIR)
        except:
            pass

        for _, filename in DATAS.items():
            out_path = os.path.join(paths.BCA_ASSAY_DIR, filename)

            if os.path.exists(out_path):
                continue

            url = DATA_URL + filename
            response = requests.get(url, stream=True)
            response.raise_for_status()

            with open(out_path, mode="wb") as f:
                for block in response.iter_content(1024):
                    f.write(block)

    def setUp(self):
        pylab.rcParams['figure.max_open_warning'] = 0

        paths.set_base_dir(os.path.abspath("."))
        self._fetch_data()

    def test_bca_assay(self):
        total_protein, _ = bca.interpret_bca_assay(
            DATAS["CK-p25"],
            samples=[
                ("3146 Cortex", 10, "A4:B6"),
                ("3131 Cortex", 10, "C4:D6"),
                ("3131 Hippocampus", 10, "E4:F6"),
                ("3131 Cerebellum", 10, "G4:H6"),
                ("3145 Cerebellum", 10, "A7:B9"),
                ("3145 Hippocampus", 10, "C7:D9"),
                ("3145 Cortex", 10, "E7:F9"),
            ],
            volumes=3000,
        )

        self.assertLess(
            abs(total_protein["3146 Cortex"] - 8000),
            500,
        )
        self.assertLess(
            abs(total_protein["3131 Cortex"] - 9000),
            500,
        )
        self.assertLess(
            abs(total_protein["3131 Hippocampus"] - 1500),
            500,
        )
        self.assertLess(
            abs(total_protein["3131 Cerebellum"] - 3500),
            500,
        )
        self.assertLess(
            abs(total_protein["3145 Cerebellum"] - 3000),
            500,
        )
        self.assertLess(
            abs(total_protein["3145 Hippocampus"] - 1500),
            500,
        )
        self.assertLess(
            abs(total_protein["3145 Cortex"] - 8500),
            500,
        )
