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
                ("base", ["low"]),
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
            psms=pd.DataFrame(columns=[
                "Proteins",
                "Sequence",
                "Modifications",
                "Validated",
                "First Scan",
                "IonScore",
                "126",
                "127",
                "128",
                "129",
                "130",
                "131",
                "Scan Paths",
            ]),
            channels=self.channels,
            groups=self.groups,
        )

        self.data.psms = self.data.psms.append(
            pd.Series(insert),
            ignore_index=True,
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

    def test_inter_normalization(self):
        psms = self.data.inter_normalize("norm")

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

    def test_inter_norm_merge(self):
        psms = self.data.inter_normalize("norm")

        cp = psms.copy()

        for chan in cp.channels.values():
            cp.psms[chan] = cp.psms[chan] * 4
            weight = "{}_weight".format(chan)
            cp.psms[weight] = cp.psms[weight] * 2

        psms = data_sets.merge_data([
            psms,
            cp,
        ])

        self.assertEqual(
            psms.psms.iloc[0]["low1"], 3,
        )
        self.assertEqual(
            psms.psms.iloc[0]["low2"], 3,
        )
        self.assertTrue(
            np.isnan(psms.psms.iloc[0]["low3"]),
        )
        self.assertEqual(
            psms.psms.iloc[0]["med"], 12,
        )
        self.assertEqual(
            psms.psms.iloc[0]["high"], 12,
        )
