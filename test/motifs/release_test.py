from test import StageTest
from test import skip_without_matlab

from pymotifs.motifs.release import Loader


class LoopsTest(StageTest):
    loader_class = Loader

    def test_can_get_best_chains_for_pdbs(self):
        val = self.loader.best_chains(['1GID', '1S72'])
        ans = {
            '1GID': set(['A', 'B']),
            '1S72': set(['0', '9'])
        }
        self.assertEquals(ans, val)

    def test_knows_all_valid_loops(self):
        val = self.loader.valid_loops('IL', ['1GID', '1S72'])
        self.assertEquals(118, len(val))

    def test_respects_the_blacklisted_loops(self):
        val = self.loader.valid_loops('HL', ['3V2F'])
        self.assertEquals(66, len(val))

    def test_can_get_all_loops_in_best_chains(self):
        val = self.loader.loops('IL', ['1GID', '1S72'])
        self.assertEquals(118, len(val))

    def test_can_get_loops_and_respect_blacklist(self):
        val = self.loader.loops('HL', ['3V2F'])
        self.assertEquals(66, len(val))

    def test_can_filter_to_just_representative_pdbs(self):
        val = self.loader.select_structures(['1GID', '2AW7', '1J5E', '1FJG'])
        ans = ['1FJG', '1GID', '2AW7']
        self.assertEquals(ans, sorted(val))


@skip_without_matlab
class ClusteringTests(StageTest):
    loader_class = Loader
