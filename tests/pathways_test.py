
from unittest import TestCase

from pyproteome.analysis import pathways


class PathwaysTest(TestCase):
    def test_gskb_pathways(self):
        for species in ["Mus musculus"]:
            gene_sets = pathways.get_gskb_pathways(species)

            for col in ["name", "set"]:
                self.assertIn(col, gene_sets.columns)

            self.assertGreater(gene_sets.shape[0], 0)

    def test_pathway_common(self):
        for species in ["Homo sapiens"]:
            gene_sets = pathways.get_pathway_common(species)

            for col in ["name", "set"]:
                self.assertIn(col, gene_sets.columns)

            self.assertGreater(gene_sets.shape[0], 0)

    def test_wikipathways(self):
        for species in ["Homo sapiens", "Mus musculus"]:
            gene_sets = pathways.get_wikipathways(species)

            for col in ["name", "set"]:
                self.assertIn(col, gene_sets.columns)

            self.assertGreater(gene_sets.shape[0], 0)

    def test_phosphomap_data(self):
        mapping = pathways.get_phosphomap_data()

        for col in ["SITE_GRP_ID", "ACC_ID", "MOD_RSD", "ORGANISM"]:
            self.assertIn(col, mapping.columns)

    def test_phosphosite(self):
        for species in ["Homo sapiens", "Mus musculus"]:
            gene_sets = pathways.get_phosphosite(species)

            for col in ["name", "set"]:
                self.assertIn(col, gene_sets.columns)

            self.assertGreater(gene_sets.shape[0], 0)

    def test_phosphosite_remap(self):
        for species in ["Homo sapiens", "Mus musculus"]:
            gene_sets = pathways.get_phosphosite_remap(species)

            for col in ["name", "set"]:
                self.assertIn(col, gene_sets.columns)

            self.assertGreater(gene_sets.shape[0], 0)

    def test_pathways(self):
        for p_sites in [True, False]:
            for species in ["Homo sapiens", "Mus musculus"]:
                gene_sets = pathways.get_pathways(species, p_sites=p_sites)

                for col in ["name", "set"]:
                    self.assertIn(col, gene_sets.columns)

                self.assertGreater(gene_sets.shape[0], 0)
