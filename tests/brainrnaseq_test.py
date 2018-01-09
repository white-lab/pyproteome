from unittest import TestCase

import brainrnaseq as brs

import pandas as pd


class BrainRNASeqTest(TestCase):
    def test_mapping_data(self):
        for species in ["Human", "Mouse"]:
            for attempt in range(3):
                map = brs.cache.get_mapping_data(
                    species=species,
                )

                self.assertIsInstance(
                    map,
                    pd.DataFrame,
                )

    def test_mapping(self):
        for attempt in range(3):
            self.assertEqual(
                brs.mapping.get_mapping(
                    gene="Jak2",
                    species="Mouse",
                ),
                "Jak2",
            )

            self.assertEqual(
                brs.mapping.get_mapping(
                    gene="JAK2",
                    species="Human",
                ),
                "JAK2",
            )

    def test_enrichment_table(self):
        for attempt in range(3):
            tab = brs.enrichments.build_enrichment_table()

            self.assertIsInstance(
                tab,
                dict,
            )

    def test_enrichments(self):
        for species in ["Human", "Mouse"]:
            for attempt in range(3):
                enrich = brs.enrichments.get_enrichments(
                    species=species,
                )

                self.assertIsInstance(
                    enrich,
                    dict,
                )
