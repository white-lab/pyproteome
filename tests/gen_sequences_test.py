from unittest import TestCase

from pycamv import gen_sequences


class GenSequencesTest(TestCase):
    def test_gen_possible_seq_test(self):
        seqs = list(
            gen_sequences.gen_possible_seq(
                "IEFTTER",
                [(1, "Phospho", ["T"])],
            )
        )
        self.assertEqual(len(seqs), 1)
        self.assertIn(
            [
                ("I", None),
                ("E", None),
                ("F", None),
                ("T", "Phospho"),
                ("T", None),
                ("E", None),
                ("R", None),
            ],
            seqs,
        )
        self.assertIn(
            [
                ("I", None),
                ("E", None),
                ("F", None),
                ("T", None),
                ("T", "Phospho"),
                ("E", None),
                ("R", None),
            ],
            seqs,
        )
