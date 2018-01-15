from unittest import TestCase

from pyproteome import modification, motif, protein, sequence


class MotifTest(TestCase):
    def setUp(self):
        self.motif = motif.Motif("O..x.-+")

    def test_repr(self):
        self.assertEqual(
            repr(self.motif),
            "<pyproteome.Motif: O..x.-+>",
        )

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
        self.sequence = sequence.Sequence(
            pep_seq="GEPNVsyICSR",
            protein_matches=[
                sequence.ProteinMatch(
                    protein=protein.Protein(
                        accession="Q9WV60",
                        gene="Gsk3b",
                        description="Glycogen synthase kinase-3 beta",
                        full_sequence=(
                            "MSGRPRTTSFAESCKPVQQPSAFGSMKVSRDKDGSKVTTVVATPGQGPD"
                            "RPQEVSYTDTKVIGNGSFGVVYQAKLCDSGELVAIKKVLQDKRFKNREL"
                            "QIMRKLDHCNIVRLRYFFYSSGEKKDEVYLNLVLDYVPETVYRVARHYS"
                            "RAKQTLPVIYVKLYMYQLFRSLAYIHSFGICHRDIKPQNLLLDPDTAVL"
                            "KLCDFGSAKQLVRGEPNVSYICSRYYRAPELIFGATDYTSSIDVWSAGC"
                            "VLAELLLGQPIFPGDSGVDQLVEIIKVLGTPTREQIREMNPNYTEFKFP"
                            "QIKAHPWTKVFRPRTPPEAIALCSRLLEYTPTARLTPLEACAHSFFDEL"
                            "RDPNVKLPNGRDTPALFNFTTQELSSNPPLATILIPPHARIQAAASPPA"
                            "NATAASDTNAGDRGQTNNAASASASNST"
                        ),
                    ),
                    rel_pos=209,
                    exact=True,
                ),
            ],
        )
        self.sequence.modifications = modification.Modifications(
            [
                # S215-p
                modification.Modification(
                    rel_pos=5,
                    mod_type="Phospho",
                    sequence=self.sequence,
                ),
                # Y216-p
                modification.Modification(
                    rel_pos=6,
                    mod_type="Phospho",
                    sequence=self.sequence,
                ),
            ],
        )

    def test_n_mers(self):
        nmers = list(
            motif.generate_n_mers(
                [self.sequence],
                letter_mod_types=[(None, "Phospho")]
            )
        )

        self.assertEqual(len(nmers), 2)
        self.assertTrue(
            all(len(i) == 15 for i in nmers)
        )
        self.assertIn(
            "VRGEPNVsYICSRYY",
            nmers,
        )
        self.assertIn(
            "RGEPNVSyICSRYYR",
            nmers,
        )

        nmers_no_filter = list(
            motif.generate_n_mers(
                [self.sequence],
            )
        )
        self.assertEqual(nmers, nmers_no_filter)


class MotifEnrichmentTest(TestCase):
    """
    Test that the motif_enrichment function runs without error on a simple
    list of sequences.
    """
    def setUp(self):
        self.sequence = sequence.Sequence(
            pep_seq="GEPNVsyICSR",
            protein_matches=[
                sequence.ProteinMatch(
                    protein=protein.Protein(
                        accession="Q9WV60",
                        gene="Gsk3b",
                        description="Glycogen synthase kinase-3 beta",
                        full_sequence=(
                            "MSGRPRTTSFAESCKPVQQPSAFGSMKVSRDKDGSKVTTVVATPGQGPD"
                            "RPQEVSYTDTKVIGNGSFGVVYQAKLCDSGELVAIKKVLQDKRFKNREL"
                            "QIMRKLDHCNIVRLRYFFYSSGEKKDEVYLNLVLDYVPETVYRVARHYS"
                            "RAKQTLPVIYVKLYMYQLFRSLAYIHSFGICHRDIKPQNLLLDPDTAVL"
                            "KLCDFGSAKQLVRGEPNVSYICSRYYRAPELIFGATDYTSSIDVWSAGC"
                            "VLAELLLGQPIFPGDSGVDQLVEIIKVLGTPTREQIREMNPNYTEFKFP"
                            "QIKAHPWTKVFRPRTPPEAIALCSRLLEYTPTARLTPLEACAHSFFDEL"
                            "RDPNVKLPNGRDTPALFNFTTQELSSNPPLATILIPPHARIQAAASPPA"
                            "NATAASDTNAGDRGQTNNAASASASNST"
                        ),
                    ),
                    rel_pos=209,
                    exact=True,
                ),
            ],
        )
        self.sequence.modifications = modification.Modifications(
            [
                # S215-p
                modification.Modification(
                    rel_pos=5,
                    mod_type="Phospho",
                    sequence=self.sequence,
                ),
                # Y216-p
                modification.Modification(
                    rel_pos=6,
                    mod_type="Phospho",
                    sequence=self.sequence,
                ),
            ],
        )
        self.sequences = list(motif.generate_n_mers([self.sequence]))
        self.foreground = self.sequences
        self.background = self.sequences

    def test_simple_enrichment(self):
        df = motif.motif_enrichment(self.foreground, self.background)
        self.assertIsNotNone(df)
