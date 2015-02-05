from nose import SkipTest

from test import StageTest

from pymotifs.species_mapping import Loader


class SpeciesMappingTest(StageTest):
    loader_class = Loader

    def test_can_generate_species_info(self):
        raise SkipTest()

    def test_generates_a_context(self):
        raise SkipTest()

    def test_downloads_a_context(self):
        raise SkipTest()

    def test_can_find_all_missing_ids(self):
        raise SkipTest()
