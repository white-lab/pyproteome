from unittest import TestCase

from pycamv import fragments


class FragmentIonsTest(TestCase):
    def test_fragment_ions(self):
        pep_seq = [
            ("N-term", "TMT6plex"),
            ("S", ()),
            ("V", ()),
            ("Y", ("Phospho",)),
            ("T", ()),
            ("E", ()),
            ("I", ()),
            ("K", ("TMT6plex")),
            ("C-term", ()),
        ]
        frag_ions = fragments.fragment_ions(
            pep_seq, 3,
        )

        hits = {
            "a_{1}^{+}": 289.21,
            "a_{2}^{+}": 388.28,
            "b_{1}^{+}": 317.20,
            "b_{2}^{+}": 416.27,
            "b_{3}^{+}": 659.30,
            "b_{4}^{+}": 760.35,
            "b_{5}^{+}": 889.39,
            "b_{6}^{+}": 1002.47,
            "y_{1}^{+}": 376.28,
            "y_{2}^{+}": 489.36,
            "y_{3}^{+}": 640.39,
            "y_{4}^{+}": 719.45,
            "y_{5}^{+}": 962.48,
            "y_{6}^{+}": 1061.55,
            "y_{7}^{+}": 1148.58,
            "pY": 216.04,
            "126": 126.13,
            "127": 127.13,
            "128": 128.13,
            "129": 129.13,
            "130": 130.14,
            "131": 131.14,
            "MH^{+2}": 689.37,
            "MH^{+2}-HPO_4": 649.39,
            "MH^{+2}-H_2O": 680.37,
            "MH^{+2}-H_2O-HPO_3": 640.39,
        }

        for name, mz in hits.items():
            self.assertIn(name, frag_ions)
            self.assertLess(abs(frag_ions[name] - mz), 0.01)
