from unittest import TestCase

import pyproteome


class MotifTest(TestCase):
    def setUp(self):
        self.motif = pyproteome.Motif("O..x.-+")

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
