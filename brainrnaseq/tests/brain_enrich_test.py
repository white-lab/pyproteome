from unittest import TestCase

import brainrnaseq


class BrainEnrichTest(TestCase):
    def test_enrich(self):
        enrich = brainrnaseq.get_enrichments("Human")
        self.assertEqual(
            enrich["INPP5D"][0], "Microglia",
        )
        self.assertEqual(
            enrich["SYT4"][0], "Neuron",
        )
        self.assertEqual(
            enrich["AGT"][0], "Astrocyte",
        )
        self.assertEqual(
            enrich["CD34"][0], "Endothelia",
        )
