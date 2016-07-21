import os

import pytest

from test import StageTest
from test import skip_without_matlab

from pymotifs import core
from pymotifs.motifs.release import Loader


class IfeTest(StageTest):
    loader_class = Loader

    def test_it_loads_representative_ifes(self):
        ifes = self.loader.ifes('1.0', ['1FJG', '1J5E', '4V4Q', '4V8I'])
        assert ifes == ['1J5E|1|A']

    def test_it_only_uses_xray_ifes(self):
        with pytest.raises(core.InvalidState):
            self.loader.ifes('1.0', ['4V6M'])

    @pytest.mark.skip("Not sure what data to use for this")
    def test_it_only_uses_structured_ifes(self):
        pass


class LoopsTest(StageTest):
    loader_class = Loader

    def loops(self, loop_type, *ifes, **kwargs):
        loop = kwargs.get('loop', '0.4')
        return self.loader.loops(loop, loop_type, ifes)

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
        loops = self.loops('HL', '1FJG|1|A', '4V4Q|1|CA')
        assert 'HL_2IL9_001' not in loops
        assert len(loops) == 64

    def test_respects_the_blacklisted_loops(self):
        loops = self.loops('HL', '2IL9|1|A', '2IL9|1|M')
        assert 'HL_2IL9_002' not in loops
        assert 'HL_2IL9_005' not in loops
        assert len(loops) == 2


class DataTests(StageTest):
    loader_class = Loader

    @skip_without_matlab
    def test_it_creates_entries_for_hl_and_il(self):
        data = self.loader.data(['1GID', '4V4Q', '1S72'], dry_run=True)
        assert len(data) == 2

    @skip_without_matlab
    def test_it_can_cluster_motifs(self):
        data = self.loader.data(['1GID', '4V4Q', '1S72'])
        assert len(data) == 2
        assert os.path.exists(data[0]['description'])
        assert os.path.exists(data[1]['description'])
