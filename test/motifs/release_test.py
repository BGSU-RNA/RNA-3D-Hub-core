import os

from test import StageTest
from test import skip_without_matlab

from pymotifs.motifs.release import Loader


class IfeTest(StageTest):
    loader_class = Loader

    def test_it_loads_representative_ifes(self):
        ifes = self.loader.ifes('1.0', ['1FJG', '1J5E', '4V4Q', '4V8I'])
        assert ifes == [
            '1FJG|1|A',
            '4V4Q|1|CA',
            '4V4Q|1|DB',
        ]

    def test_it_only_uses_xray_ifes(self):
        ifes = self.loader.ifes('1.0', ['4V6M'])
        assert ifes == []

    def test_it_only_uses_structured_ifes(self):
        ifes = self.loader.ifes('1.0', ['1FJG'])
        assert '1FJG|1|X' not in ifes
        assert ifes == ['1FJG|1|A']


class LoopsTest(StageTest):
    loader_class = Loader

    def loops(self, loop_type, *ifes):
        return self.loader.loops('0.1', loop_type, ifes)

    def test_uses_all_valid_loops(self):
        assert len(self.loops('IL', '1FJG|1|A')) == 62
        assert len(self.loops('HL', '1FJG|1|A')) == 32

    def test_it_uses_loops_that_passed_loop_qa(self):
        loops = self.loops('HL', '1S72|1|0', '1S72|1|9')
        assert 'HL_1S72_004' not in loops
        assert 'HL_1S72_018' not in loops
        assert 'HL_1S72_021' not in loops
        assert len(loops) == 63

    def test_it_uses_loops_from_ife_chains(self):
        loops = self.loop('HL', '1FJG|1|A', '4V4Q|1|CA')
        assert '' not in loops
        assert len(loops) == 10000

    def test_respects_the_blacklisted_loops(self):
        loops = self.loops('HL', '2IL9')
        assert 'HL_2IL9_002' not in loops
        assert 'HL_2IL9_005' not in loops
        assert len(loops) == 66


class DataTests(StageTest):
    loader_class = Loader

    def test_it_creates_entries_for_hl_and_il(self):
        data = self.loader.data(['1GID', '4V4Q', '1S72'], dry_run=True)
        assert len(data) == 2

    @skip_without_matlab
    def test_it_can_cluster_motifs(self):
        data = self.loader.data(['1GID', '4V4Q', '1S72'])
        assert len(data) == 2
        assert os.path.exists(data[0]['description'])
        assert os.path.exists(data[1]['description'])
