from unittest import TestCase

import brainrnaseq as brs

import pandas as pd


class BrainRNASeqTest(TestCase):
    def test_mapping_data(self):
        items = {
            "Homo sapiens": ["A1BG", "A1B|ABG|GAB|HYST2477", "16S rRNA", "-"],
            "Mus musculus": ["Pzp", "A1m|A2m|AI893533|MAM", "ND2", "-"],
        }

        for species, vals in items.items():
            for _ in range(2):
                map = brs.cache.get_mapping_data(
                    species=species,
                )

                self.assertIsInstance(
                    map,
                    pd.DataFrame,
                )
                self.assertEqual(
                    map.index[0],
                    vals[0],
                )
                self.assertEqual(
                    map.iloc[0]["Synonyms"],
                    vals[1],
                )
                self.assertEqual(
                    map.index[-1],
                    vals[2],
                )
                self.assertEqual(
                    map.iloc[-1]["Synonyms"],
                    vals[3],
                )

            map = brs.cache.get_mapping_data(
                species=species,
                force=True,
            )
            self.assertIsInstance(
                map,
                pd.DataFrame,
            )

    def test_mapping(self):
        for _ in range(2):
            for gene in ["Jak2", "Fd17"]:
                self.assertEqual(
                    brs.mapping.get_symbol_mapping(
                        gene=gene,
                        species="Mus musculus",
                    ),
                    "Jak2",
                )

            self.assertEqual(
                brs.mapping.get_symbol_mapping(
                    gene="NaNNaNNaN",
                    species="Mus musculus",
                ),
                None,
            )

            for gene in ["JAK2", "JTK10", "THCYT3"]:
                self.assertEqual(
                    brs.mapping.get_symbol_mapping(
                        gene=gene,
                        species="Homo sapiens",
                    ),
                    "JAK2",
                )

            # Ambiguous gene
            self.assertEqual(
                brs.mapping.get_symbol_mapping(
                    gene="AP-1",
                    species="Homo sapiens",
                ),
                "FOS",
            )

    def test_barres_seq_data(self):
        brs.cache.get_barres_seq_data()

        items = {
            "Homo sapiens": ["Gene", "1/2-SBSRNA4", "ZZZ3"],
            "Mus musculus": ["gene", "0610005C13Rik", "Zzz3"],
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
        for _ in range(2):
            tab = brs.enrichments.build_enrichment_table()

            self.assertIsInstance(
                tab,
                dict,
            )

    def test_enrichments(self):
        items = {
            "Homo sapiens": {
                "AGT": "Astrocyte",
                "GFAP": "Astrocyte",
                "IDI2-AS1": "Astrocyte",
                "CD34": "Endothelia",
                "RNU11": "Endothelia",
                "CCL3L1": "Microglia",
                "CD33": "Microglia",
                "INPP5D": "Microglia",
                "SIGLEC5": "Microglia",
                "SIGLEC8": "Microglia",
                "SIGLEC9": "Microglia",
                "SIGLEC10": "Microglia",
                "SIGLEC14": "Microglia",
                "FOLH1": "Myelinating Oligodendrocytes",
                "CD22": "Myelinating Oligodendrocytes",
                "MAG": "Myelinating Oligodendrocytes",
                "GAD2": "Neuron",
                "SYT4": "Neuron",
            },
            "Mus musculus": {
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
        for species, vals in items.items():
            for _ in range(2):
                enrich = brs.enrichments.get_enrichments(
                    species=species,
                )

                self.assertIsInstance(
                    enrich,
                    dict,
                )
                for key, val in vals.items():
                    self.assertEqual(
                        enrich[key],
                        val,
                    )
