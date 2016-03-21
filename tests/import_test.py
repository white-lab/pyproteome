from unittest import TestCase


class ImportTest(TestCase):
    def test_imports(self):
        from pycamv import (
            compare, fragments, gen_sequences, mascot, masses, ms_labels,
            proteowizard, scans, validate,
        )
        from pyproteome import (
            analysis, bca, camv, data_sets, enrichment, fetch_data, levels,
            loading, logo, modification, motif, paths, protein, sequence,
            utils, version,
        )
