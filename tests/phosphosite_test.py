
from unittest import TestCase

# import pylab

from pyproteome import phosphosite


class PhosphositeTest(TestCase):
    # def setUp(self):
    #     pylab.rcParams['figure.max_open_warning'] = 0

    def test_generate_logos(self):
        for species in ["Human", "Mouse"]:
            phosphosite.generate_logos(
                species,
                min_foreground=50,
            )
