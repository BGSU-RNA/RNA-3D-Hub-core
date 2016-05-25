import pytest

import collections as coll

from test import StageTest

from pymotifs.interactions.summary import Loader


class QueryTest(StageTest):
    loader_class = Loader

    def test_it_knows_if_no_data_exists(self):
        assert self.loader.has_data('0GID') is False

    @pytest.mark.xfail(reason='Must bootstrap')
    def test_knows_if_data_exists(self):
        assert self.loader.has_data('0GID') is True


class BasicCountsTest(StageTest):
    loader_class = Loader

    def test_it_can_increment_bps(self):
        data = coll.defaultdict(int)
        val = dict(self.loader.increment_bp(data, 'cWW', 3))
        assert val == {
            'cWW': 1,
            'bps': 1,
            'total': 1,
            'lr_cWW': 0,
            'lr_total': 0,
            'lr_bps': 0,
        }

    def test_it_can_increment_long_range_bps(self):
        data = coll.defaultdict(int)
        val = dict(self.loader.increment_bp(data, 'cWW', 5))
        assert val == {
            'cWW': 1,
            'bps': 1,
            'total': 1,
            'lr_cWW': 1,
            'lr_total': 1,
            'lr_bps': 1,
        }

    def test_it_can_increment_stacks(self):
        data = coll.defaultdict(int)
        val = dict(self.loader.increment_stacks(data, 's33', 3))
        assert val == {
            's33': 1,
            'stacks': 1,
            'total': 1,
            'lr_s33': 0,
            'lr_total': 0,
            'lr_stacks': 0,
        }

    def test_it_can_increment_long_range_stacks(self):
        data = coll.defaultdict(int)
        val = dict(self.loader.increment_stacks(data, 's33', 5))
        assert val == {
            's33': 1,
            'stacks': 1,
            'total': 1,
            'lr_s33': 1,
            'lr_total': 1,
            'lr_stacks': 1,
        }

    def test_it_can_increment_bphs(self):
        data = coll.defaultdict(int)
        val = dict(self.loader.increment_bphs(data, 1, 2, '1BPh', 3))
        assert val == {
            '1BPh': 1,
            'bphs': 1,
            'total': 1,
            'lr_1BPh': 0,
            'lr_total': 0,
            'lr_bphs': 0,
        }

    def test_it_can_increment_long_range_bphs(self):
        data = coll.defaultdict(int)
        val = dict(self.loader.increment_bphs(data, 1, 2, '1BPh', 5))
        assert val == {
            '1BPh': 1,
            'bphs': 1,
            'total': 1,
            'lr_1BPh': 1,
            'lr_total': 1,
            'lr_bphs': 1,
        }

    def test_it_does_not_count_self_0BPh(self):
        data = coll.defaultdict(int)
        val = dict(self.loader.increment_bphs(data, 1, 1, '1BPh', 5))
        assert val == {}

    def test_it_does_not_count_near_bp(self):
        data = coll.defaultdict(int)
        val = dict(self.loader.increment_bp(data, 'ncWW', 5))
        assert val == {}

    def test_it_does_not_count_near_bphs(self):
        data = coll.defaultdict(int)
        val = dict(self.loader.increment_bphs(data, 1, 2, 'n1BPh', 5))
        assert val == {}

    def test_it_does_not_count_near_stacks(self):
        data = coll.defaultdict(int)
        val = dict(self.loader.increment_stacks(data, 'ns55', 5))
        assert val == {}

    def test_it_does_not_increment_with_no_bp(self):
        data = coll.defaultdict(int)
        val = dict(self.loader.increment_bp(data, None, 5))
        assert val == {}

    def test_it_does_not_increment_with_no_stack(self):
        data = coll.defaultdict(int)
        val = dict(self.loader.increment_stacks(data, None, 5))
        assert val == {}

    def test_it_does_not_increment_with_no_bphs(self):
        data = coll.defaultdict(int)
        val = dict(self.loader.increment_bphs(data, 1, 1, None, 5))
        assert val == {}


class DataTest(StageTest):
    loader_class = Loader

    def data(self, pdb):
        data = {}
        for entry in self.loader.data(pdb):
            data[entry['unit_id']] = entry
        return data

    def test_it_counts_all_nts_in_the_structure(self):
        val = self.data('4V4Q')
        assert len(val) == 8976

    def test_it_only_counts_for_if_unit_is_first(self):
        val = self.data('4V4Q')
        assert val['4V4Q|1|BA|A|29']['total'] == 4

    def test_it_computes_the_correct_counts(self):
        val = self.data('4V4Q')
        assert val['4V4Q|1|DB|C|645'] == {
            'unit_id': '4V4Q|1|DB|C|645',
            'pdb_id': '4V4Q',
            'model': 1,
            'chain': 'DB'
        }
        assert val['4V4Q|1|BA|A|29']['total'] == 4
