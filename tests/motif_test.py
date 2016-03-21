from unittest import TestCase

import pyproteome as pyp


class MotifTest(TestCase):
    def setUp(self):
        self.motif = pyp.Motif("O..x.-+")

    def test_match(self):
        self.assertIn("IEFtFER", self.motif)
        self.assertIn("LEFsFER", self.motif)
        self.assertIn("VEFtFER", self.motif)
        self.assertIn("MEFsFER", self.motif)

        self.assertIn("IEFsFER", self.motif)
        self.assertIn("IEFtFER", self.motif)

        self.assertIn("IEFsFER", self.motif)
        self.assertIn("IEFsFEK", self.motif)

    def test_no_match(self):
        self.assertNotIn("QEFtFER", self.motif)
        self.assertNotIn("HEFsFER", self.motif)
        self.assertNotIn("EEFtFER", self.motif)
        self.assertNotIn("DEFsFER", self.motif)

        self.assertNotIn("IEFYFER", self.motif)
        self.assertNotIn("IEFSFER", self.motif)
        self.assertNotIn("IEFTFER", self.motif)

        self.assertNotIn("IEFsFED", self.motif)
        self.assertNotIn("IEFsFEE", self.motif)

    def test_match_motif(self):
        # Motifs will match themselves and any more specific motifs, though
        # they will not match less specific motifs
        self.assertTrue(self.motif.match("O..x.-+"))
        self.assertTrue(self.motif.match("O..x.E+"))
        self.assertTrue(self.motif.match("O..s.-+"))
        self.assertTrue(self.motif.match("O..t.-+"))
        self.assertTrue(self.motif.match("M..x.-+"))
        self.assertTrue(self.motif.match("O..x.-K"))

        self.assertFalse(self.motif.match("...x.-+"))
        self.assertFalse(self.motif.match("O....-+"))
        self.assertFalse(self.motif.match("O..x..+"))
        self.assertFalse(self.motif.match("O..x.-."))


class GenerateNMersTest(TestCase):
    def setUp(self):
        self.sequence = pyp.Sequence(
            pep_seq="GEPNVsyICSR",
            protein_matches=[
                pyp.ProteinMatch(
                    protein=pyp.Protein(
                        accession="Q9WV60",
                    ),
                    rel_pos=209,
                    exact=True,
                ),
            ],
        )
        self.sequence.modifications = pyp.Modifications(
            [
                # S215-p
                pyp.Modification(
                    rel_pos=5,
                    mod_type="Phospho",
                    sequence=self.sequence,
                ),
                # Y216-p
                pyp.Modification(
                    rel_pos=6,
                    mod_type="Phospho",
                    sequence=self.sequence,
                ),
            ],
        )

    def test_n_mers(self):
        nmers = list(
            pyp.motif.generate_n_mers(
                [self.sequence],
                letter_mod_types=[(None, "Phospho")]
            )
        )

        self.assertEqual(len(nmers), 2)
        self.assertTrue(
            all(len(i) == 15 for i in nmers)
        )
        self.assertEqual(
            nmers[0],
            "VRGEPNVsYICSRYY",
        )
        self.assertEqual(
            nmers[1],
            "RGEPNVSyICSRYYR",
        )

        nmers_no_filter = list(
            pyp.motif.generate_n_mers(
                [self.sequence],
            )
        )
        self.assertEqual(nmers, nmers_no_filter)


class MotifEnrichmentTest(TestCase):
    """
    Test that the motif_enrichment function runs without error on a simple
    set of Sequence objects.
    """
    def setUp(self):
        self.sequence = pyp.Sequence(
            pep_seq="GEPNVsyICSR",
            protein_matches=[
                pyp.ProteinMatch(
                    protein=pyp.Protein(
                        accession="Q9WV60",
                    ),
                    rel_pos=209,
                    exact=True,
                ),
            ],
        )
        self.sequence.modifications = pyp.Modifications(
            [
                # S215-p
                pyp.Modification(
                    rel_pos=5,
                    mod_type="Phospho",
                    sequence=self.sequence,
                ),
                # Y216-p
                pyp.Modification(
                    rel_pos=6,
                    mod_type="Phospho",
                    sequence=self.sequence,
                ),
            ],
        )
        self.sequences = [self.sequence]
        self.foreground = self.sequences
        self.background = self.sequences

    def test_simple_enrichment(self):
        df = pyp.motif.motif_enrichment(self.foreground, self.background)
        self.assertIsNotNone(df)
