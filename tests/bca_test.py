
import os
from unittest import TestCase

import pylab

from pyproteome import bca, paths

from . import utils


DATA_URL = "https://github.com/white-lab/pyproteome-data/raw/master/test_data/"

DATAS = {
    "CK-p25": "CK-p25-bca.xlsx",
}


class BCATest(TestCase):
    def setUp(self):
        paths.set_base_dir(os.path.abspath("."))

        utils.fetch_data(
            dir=paths.BCA_ASSAY_DIR,
            datas=DATAS,
            base_url=DATA_URL,
        )

        pylab.rcParams['figure.max_open_warning'] = 0

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
            abs(total_protein["3146 Cortex"][0] - 8000),
            500,
        )
        self.assertLess(
            abs(total_protein["3131 Cortex"][0] - 9000),
            500,
        )
        self.assertLess(
            abs(total_protein["3131 Hippocampus"][0] - 1500),
            500,
        )
        self.assertLess(
            abs(total_protein["3131 Cerebellum"][0] - 3500),
            500,
        )
        self.assertLess(
            abs(total_protein["3145 Cerebellum"][0] - 3000),
            500,
        )
        self.assertLess(
            abs(total_protein["3145 Hippocampus"][0] - 1500),
            500,
        )
        self.assertLess(
            abs(total_protein["3145 Cortex"][0] - 8500),
            500,
        )
