from unittest import TestCase


class ImportTest(TestCase):
    def test_imports(self):
        from pyproteome import (
            analysis, bca, camv, data_sets, enrichment, fetch_data, levels,
            loading, logo, modification, motif, paths, protein, sequence,
            utils, version,
        )
