
from __future__ import division

from unittest import TestCase

from pyproteome import data_sets, loading, modification, protein

import pandas as pd
import numpy as np
from collections import OrderedDict


class NormalizationTest(TestCase):
    def setUp(self):
        self.prots = protein.Proteins(
            proteins=[
                protein.Protein(
                    accession="P03995"  # GFAP_MOUSE
                )
            ],
        )
        self.seq = loading.extract_sequence(self.prots, "QEADEATLAR")

        self.mods = modification.Modifications(
            mods=[
                modification.Modification(
                    rel_pos=0,
                    mod_type="TMT6plex",
                    nterm=True,
                    sequence=self.seq,
                ),
            ],
        )

        self.seq.modifications = self.mods
        self.channels = OrderedDict(
            [
                ("low1", "126"),
                ("low2", "127"),
                ("low3", "128"),
                ("med", "129"),
                ("high", "130"),
                ("norm", "131"),
            ]
        )
        self.groups = OrderedDict(
            [
                ("base", ["low1", "low2", "low3"]),
                ("stim", ["med", "high"]),
            ]
        )

        insert = {
            "Proteins": self.prots,
            "Sequence": self.seq,
            "Modifications": self.mods,
            "126": 1e4,
            "127": 1e4,
            "128": np.nan,
            "129": 4e4,
            "130": 4e4,
            "131": 1e4,
            "Validated": True,
            "First Scan": {17015},
            "Scan Paths": {"blank"},
            "IonScore": 100,
        }

        self.data = data_sets.DataSet(
            channels=self.channels,
            groups=self.groups,
        )

        self.data.psms = self.data.psms.append(
            pd.Series(insert),
            ignore_index=True,
        )

    def test_dropna(self):
        psms = self.data.dropna()
        self.assertEqual(
            len(psms.psms), 0,
        )

    def test_normalization(self):
        norm_levels = OrderedDict([
            ("126", 1),
            ("127", np.nan),
            ("128", 1),
            ("129", 2),
            ("130", 1),
            ("131", 1),
        ])
        psms = self.data.normalize(norm_levels)

        self.assertEqual(
            psms.psms.iloc[0]["126_norm"], 1e4,
        )
        self.assertTrue(
            np.isnan(psms.psms.iloc[0]["127_norm"]),
        )
        self.assertTrue(
            np.isnan(psms.psms.iloc[0]["128_norm"]),
        )
        self.assertEqual(
            psms.psms.iloc[0]["129_norm"], 2e4,
        )
        self.assertEqual(
            psms.psms.iloc[0]["130_norm"], 4e4,
        )
        self.assertEqual(
            psms.psms.iloc[0]["131_norm"], 1e4,
        )
        self.assertEqual(
            psms.psms.iloc[0]["Fold Change"], 1/3,
        )
        self.assertTrue(
            np.isnan(psms.psms.iloc[0]["p-value"]),
        )

    def test_inter_normalization(self):
        psms = self.data.inter_normalize(["norm"])

        self.assertEqual(
            psms.psms.iloc[0]["low1"], 1,
        )
        self.assertEqual(
            psms.psms.iloc[0]["low2"], 1,
        )
        self.assertTrue(
            np.isnan(psms.psms.iloc[0]["low3"]),
        )
        self.assertEqual(
            psms.psms.iloc[0]["med"], 4,
        )
        self.assertEqual(
            psms.psms.iloc[0]["high"], 4,
        )
        self.assertEqual(
            psms.psms.iloc[0]["Fold Change"], 1/4,
        )
        self.assertTrue(
            np.isnan(psms.psms.iloc[0]["p-value"]),
        )

    def test_inter_norm_merge(self):
        psms = self.data.copy()

        cp = psms.copy()

        for name, chan in cp.channels.items():
            if name == "norm":
                # cp.psms[chan] /= 2
                continue

            cp.psms[chan] *= 4
            weight = "{}_weight".format(chan)
            cp.psms[weight] = cp.psms[chan] / 2

        cp.psms["131_C"] = cp.psms[cp.channels["high"]] * 2
        cp.channels["newhigh"] = "131_C"

        psms = data_sets.merge_data(
            [
                psms,
                cp,
            ],
            norm_channels=["norm"],
        )

        self.assertEqual(
            psms.psms.iloc[0]["low1"] /
            psms.psms.iloc[0]["norm"], 3,
        )
        self.assertEqual(
            psms.psms.iloc[0]["low2"] /
            psms.psms.iloc[0]["norm"], 3,
        )
        self.assertTrue(
            np.isnan(psms.psms.iloc[0]["low3"]),
        )
        self.assertEqual(
            psms.psms.iloc[0]["med"] /
            psms.psms.iloc[0]["norm"], 12,
        )
        self.assertEqual(
            psms.psms.iloc[0]["high"] /
            psms.psms.iloc[0]["norm"], 12,
        )
        self.assertEqual(
            psms.psms.iloc[0]["newhigh"] /
            psms.psms.iloc[0]["norm"], 32,
        )
        self.assertEqual(
            psms.psms.iloc[0]["Fold Change"], 1/4,
        )
        self.assertTrue(
            np.isnan(psms.psms.iloc[0]["p-value"]),
        )

    def test_nan_norm(self):
        self.data.psms["131"] = np.nan
        psms = self.data.inter_normalize(["norm"])

        self.assertTrue(
            np.isnan(psms.psms.iloc[0]["low1"]),
        )
        self.assertTrue(
            np.isnan(psms.psms.iloc[0]["low2"]),
        )
        self.assertTrue(
            np.isnan(psms.psms.iloc[0]["low3"]),
        )
        self.assertTrue(
            np.isnan(psms.psms.iloc[0]["med"]),
        )
        self.assertTrue(
            np.isnan(psms.psms.iloc[0]["high"]),
        )
        self.assertTrue(
            np.isnan(psms.psms.iloc[0]["Fold Change"]),
        )
        self.assertTrue(
            np.isnan(psms.psms.iloc[0]["p-value"]),
        )
