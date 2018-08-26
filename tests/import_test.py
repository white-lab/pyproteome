from unittest import TestCase


class ImportTest(TestCase):
    def test_imports(self):
        import pyproteome
        import brainrnaseq

        # Fix linter warnings
        pyproteome.DUMMY = None
        brainrnaseq.DUMMY = None
