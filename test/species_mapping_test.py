from unittest import TestCase

from test import StageTest

from pymotifs.species_mapping import Parser
from pymotifs.species_mapping import Loader


class ParsingTest(TestCase):

    def parse(self, filename):
        parser = Parser()
        with open(filename, 'rb') as raw:
            return parser.parse(raw.read())

    def test_can_get_correct_data_given_species(self):
        val = self.parse('test/files/species-mapping/species.xml')
        ans = {
            'species_mapping_id': 562,
            'species_id': 562,
            'species_name': 'Escherichia coli'
        }
        assert val == ans

    def test_can_get_correct_data_given_subspecies(self):
        val = self.parse('test/files/species-mapping/subspecies.xml')
        ans = {
            'species_mapping_id': 656435,
            'species_id': 562,
            'species_name': 'Escherichia coli'
        }
        assert val == ans

    def test_can_get_correct_data_given_genus(self):
        val = self.parse('test/files/species-mapping/genus.xml')
        ans = {
            'species_mapping_id': 4930,
            'species_id': None,
            'species_name': None
        }
        assert val == ans


class SpeciesMappingTest(StageTest):
    loader_class = Loader

    def test_can_transform_to_taxon_ids(self):
        val = self.loader.to_process(['1FJ0', '1S72', '4V4Q'])
        ans = [2238, 562]
        self.assertEquals(ans, val)

    def test_knows_if_has_data(self):
        self.assertTrue(self.loader.has_data(562))

    def test_knows_if_has_no_data(self):
        self.assertFalse(self.loader.has_data(-1))

    def test_can_fetch_for_one_species(self):
        val = self.loader.data(1590)
        self.assertEquals(1590, val.species_id)
        self.assertEquals('Lactobacillus plantarum', val.species_name)
        self.assertEquals(1590, val.species_mapping_id)
