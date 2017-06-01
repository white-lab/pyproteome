from unittest import TestCase


class ImportTest(TestCase):
    def test_imports(self):
        from pyproteome import (
            analysis, bca, camv, data_sets, discoverer, enrichment, fetch_data,
            fetch_data, levels, loading, logo, modification, motif, paths,
            pride, protein, sequence, utils, version,
        )
