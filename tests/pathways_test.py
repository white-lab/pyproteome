
from unittest import TestCase

from pyproteome import pathways


class PathwaysTest(TestCase):
    def test_gskb_pathways(self):
        for species in ['Mus musculus']:
            gene_sets = pathways.gskb.get_gskb_pathways(
                species,
            )

            for col in ['name', 'set']:
                self.assertIn(col, gene_sets.columns)

            self.assertGreater(gene_sets.shape[0], 0)

    def test_pathway_common(self):
        for species in ['Homo sapiens']:
            gene_sets = pathways.pathwayscommon.get_pathway_common(
                species,
            )

            for col in ['name', 'set']:
                self.assertIn(col, gene_sets.columns)

            self.assertGreater(gene_sets.shape[0], 0)

    def test_wikipathways(self):
        for species in ['Homo sapiens', 'Mus musculus']:
            gene_sets = pathways.wikipathways.get_wikipathways(
                species,
            )

            for col in ['name', 'set']:
                self.assertIn(col, gene_sets.columns)

            self.assertGreater(gene_sets.shape[0], 0)

    def test_phosphomap_data(self):
        mapping = pathways.psp.get_phosphomap_data()

        for col in ['SITE_GRP_ID', 'ACC_ID', 'MOD_RSD', 'ORGANISM']:
            self.assertIn(col, mapping.columns)

    def test_phosphosite(self):
        for species in ['Homo sapiens', 'Mus musculus']:
            for remap in [False, True]:
                gene_sets = pathways.psp.get_phosphosite(
                    species,
                    remap=remap,
                )

                for col in ['name', 'up_set']:
                    self.assertIn(col, gene_sets.columns)

                self.assertGreater(gene_sets.shape[0], 0)

    def test_phosphosite_regulation(self):
        for species in ['Homo sapiens', 'Mus musculus']:
            for remap in [False, True]:
                gene_sets = pathways.psp.get_phosphosite_regulation(
                    species,
                    remap=remap,
                )

                for col in ['name', 'up_set', 'down_set']:
                    self.assertIn(col, gene_sets.columns)

                self.assertGreater(gene_sets.shape[0], 0)

    def test_pathways(self):
        for p_sites in [True, False]:
            for species in ['Homo sapiens', 'Mus musculus']:
                for remap in [False, True]:
                    if not p_sites and remap:
                        continue

                    gene_sets = pathways.get_pathways(
                        species,
                        p_sites=p_sites,
                        remap=remap,
                    )

                    for col in ['name'] + (
                        ['up_set', 'down_set'] if p_sites else ['set']
                    ):
                        self.assertIn(col, gene_sets.columns)

                    self.assertGreater(gene_sets.shape[0], 0)

    def test_PrPDF(self):
        pdf = pathways.enrichments.PrPDF(range(1, 11))
        self.assertEqual(
            pdf.cdf(1),
            0,
        )
        self.assertEqual(
            pdf.cdf(2),
            .1,
        )
        self.assertEqual(
            pdf.cdf(10),
            .9,
        )
        self.assertEqual(
            pdf.sf(2),
            .8,
        )
        self.assertEqual(
            pdf.sf(10),
            0,
        )
        self.assertEqual(
            pdf.pdf(2),
            .1,
        )
        self.assertEqual(
            pdf.cdf(2) + pdf.pdf(2) + pdf.sf(2),
            1,
        )
