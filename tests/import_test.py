from unittest import TestCase


class ImportTest(TestCase):
    def test_imports(self):
        from pyproteome import (
            analysis, bca, camv, data_sets, discoverer, enrichment, fetch_data,
            fetch_data, icelogo, levels, loading, logo, modification, motif,
            paths, plogo, pride, protein, sequence, utils, version, weblogo,
        )
