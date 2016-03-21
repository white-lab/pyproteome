from unittest import TestCase

from pycamv import masses


class MassesTest(TestCase):
    def test_exact_mass(self):
        self.assertAlmostEqual(
            masses.exact_mass({"H": 1}),
            1.007940721855,
        )

        self.assertAlmostEqual(
            masses.exact_mass({"H": [1]}),
            1.007825,
        )

        self.assertAlmostEqual(
            masses.exact_mass({"H": [1, 1]}),
            3.021927,
        )

        self.assertAlmostEqual(
            masses.exact_mass({"H": [0, 1]}),
            2.014102,
        )
