from unittest import TestCase

import brainrnaseq as brs

import pandas as pd


class BrainRNASeqTest(TestCase):
    def test_mapping_data(self):
        items = {
            "Human": ["A1BG", "A1B|ABG|GAB|HYST2477", "16S rRNA", "-"],
            "Mouse": ["Pzp", "A1m|A2m|AI893533|MAM", "ND2", "-"],
        }

        for species in ["Human", "Mouse"]:
            for attempt in range(2):
                map = brs.cache.get_mapping_data(
                    species=species,
                )

                self.assertIsInstance(
                    map,
                    pd.DataFrame,
                )
                self.assertEqual(
                    map.index[0],
                    items[species][0],
                )
                self.assertEqual(
                    map.iloc[0]["Synonyms"],
                    items[species][1],
                )
                self.assertEqual(
                    map.index[-1],
                    items[species][2],
                )
                self.assertEqual(
                    map.iloc[-1]["Synonyms"],
                    items[species][3],
                )

    def test_mapping(self):
        for attempt in range(2):
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

    def test_barres_seq_data(self):
        brs.cache.get_barres_seq_data()

        items = {
            "Human": ["Gene", "1/2-SBSRNA4", "ZZZ3"],
            "Mouse": ["gene", "0610005C13Rik", "Zzz3"],
        }

        for species, (col, first, last) in items.items():
            self.assertEqual(
                brs.cache.SPECIES_DATA[species].iloc[0][col],
                first,
            )
            self.assertEqual(
                brs.cache.SPECIES_DATA[species].iloc[-1][col],
                last,
            )

    def test_enrichment_table(self):
        for attempt in range(2):
            tab = brs.enrichments.build_enrichment_table()

            self.assertIsInstance(
                tab,
                dict,
            )

    def test_enrichments(self):
        items = {
            "Human": {
                "GFAP": "Astrocyte",
                "IDI2-AS1": "Astrocyte",
                "RNU11": "Endothelia",
                "CCL3L1": "Microglia",
                "CD33": "Microglia",
                "SIGLEC5": "Microglia",
                "SIGLEC8": "Microglia",
                "SIGLEC9": "Microglia",
                "SIGLEC10": "Microglia",
                "SIGLEC14": "Microglia",
                "FOLH1": "Myelinating Oligodendrocytes",
                "CD22": "Myelinating Oligodendrocytes",
                "MAG": "Myelinating Oligodendrocytes",
                "GAD2": "Neuron",
            },
            "Mouse": {
                "Sumo2": "Astrocyte",
                "AU021092": "Endothelia",
                "OncoM": "Microglia",
                "Siglec-1": "Microglia",
                "Siglec1": "Microglia",
                "Cd22": "Microglia",
                "Cd33": "Microglia",
                "Siglec-3": "Microglia",
                "mSiglec-E": "Microglia",
                "Siglece": "Microglia",
                "Siglec5": "Microglia",
                "Siglec9": "Microglia",
                "Siglech": "Microglia",
                "Siglec-H": "Microglia",
                "Siglecl1": "Microglia",
                "Siglec12": "Microglia",
                "siglec-4a": "Myelinating Oligodendrocytes",
                "Otm": "Myelinating Oligodendrocytes",
                "Reln": "Neuron",
            },
        }
        for species in ["Human", "Mouse"]:
            for attempt in range(2):
                enrich = brs.enrichments.get_enrichments(
                    species=species,
                )

                self.assertIsInstance(
                    enrich,
                    dict,
                )
                for key, val in items[species].items():
                    self.assertEqual(
                        enrich[key][0],
                        val,
                    )
