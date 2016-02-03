from unittest import TestCase

import pyproteome as pyp


class MotifTest(TestCase):
    def setUp(self):
        self.motif = pyp.Motif("O..x.-+")

    def test_match(self):
        self.assertTrue("IEFyFER" in self.motif)
        self.assertTrue("LEFyFER" in self.motif)
        self.assertTrue("VEFyFER" in self.motif)
        self.assertTrue("MEFyFER" in self.motif)

        self.assertTrue("IEFyFER" in self.motif)
        self.assertTrue("IEFsFER" in self.motif)
        self.assertTrue("IEFtFER" in self.motif)

        self.assertTrue("IEFyFER" in self.motif)
        self.assertTrue("IEFyFEK" in self.motif)

    def test_no_match(self):
        self.assertFalse("QEFyFER" in self.motif)
        self.assertFalse("HEFyFER" in self.motif)
        self.assertFalse("EEFyFER" in self.motif)
        self.assertFalse("DEFyFER" in self.motif)

        self.assertFalse("IEFYFER" in self.motif)
        self.assertFalse("IEFSFER" in self.motif)
        self.assertFalse("IEFTFER" in self.motif)

        self.assertFalse("IEFyFED" in self.motif)
        self.assertFalse("IEFyFEE" in self.motif)


class GenerateNMersTest(TestCase):
    def setUp(self):
        self.sequence = pyp.Sequence(
            pep_seq="GEPNVSyICSR",
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
