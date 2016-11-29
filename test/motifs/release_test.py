import os

import pytest

from test import StageTest
from test import skip_without_matlab

from pymotifs.motifs.release import Loader


class IfeTest(StageTest):
    loader_class = Loader

    def test_it_loads_representative_ifes(self):
        assert self.loader.ifes('1.0') == [
            '157D|1|A+157D|1|B', '1DUH|1|A', '1EIY|1|C', '1ET4|1|E',
            '1G59|1|D', '1GID|1|B', '1IBK|1|X', '1J5E|1|A', '1KOG|1|P',
            '1MDG|1|A', '1UTD|1|0', '1VY4|1|AY', '1VY4|1|BA', '1VY4|1|BB',
            '1WMQ|1|D', '1X8W|1|B', '1X8W|1|C', '2HOJ|1|A', '2IL9|1|A',
            '3CW5|1|A', '4A3G|1|P', '4CS1|1|A', '4FTE|1|R', '4PMI|1|A',
            '4Q0B|1|T', '4Q0B|1|t', '4V88|1|A5+4V88|1|A8', '4V88|1|A6',
            '4V88|1|A7', '4V9K|1|CW'
        ]

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

    @pytest.mark.skip("Do not want to run full clustering yet")
    def test_it_creates_entries_for_hl_and_il(self):
        data = self.loader.data(['1GID', '4V4Q', '1S72'], dry_run=True)
        assert len(data) == 2

    @pytest.mark.skip("Do not want to run full clustering yet")
    def test_it_can_cluster_motifs(self):
        data = self.loader.data(['1GID', '4V4Q', '1S72'])
        assert len(data) == 2
        assert os.path.exists(data[0]['description'])
        assert os.path.exists(data[1]['description'])
