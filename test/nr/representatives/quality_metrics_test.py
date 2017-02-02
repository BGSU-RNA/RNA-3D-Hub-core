from test import StageTest

from pymotifs.nr.groups.simplified import Grouper
from pymotifs.nr.representatives import QualityMetrics


class LoadingDataTests(StageTest):
    loader_class = QualityMetrics

    def test_it_can_load_quality_data_when_missing(self):
        assert self.loader.load_quality([{'pdb': '0S72'}]) == [{
            'pdb': '0S72',
            'quality': {
                'has': set(),
                'rsrz': 100,
                'backbone': 100,
                'clashscore': 500,
            }
        }]

    def test_it_can_load_known(self):
        assert self.loader.load_quality([{'pdb': '1FJG'}]) == [{
            'pdb': '1FJG',
            'quality': {
                'has': set(['rsrz', 'clashscore', 'backbone']),
                'rsrz': 5.82,
                'backbone': 7.0,
                'clashscore': 29.48,
            }
        }]

    def test_it_can_load_when_missing_partial(self):
        assert self.loader.load_quality([{'pdb': '1GID'}]) == [{
            'pdb': '1GID',
            'quality': {
                'has': set([u'clashscore']),
                'rsrz': 100,
                'backbone': 100,
                'clashscore': 36.87,
            }
        }]

    def test_can_load_all_data(self):
        members = [
            {'pdb': '1GID'},
            {'pdb': '0S72'},
            {'pdb': '1FJG'},
        ]
        assert self.loader.load_quality(members) == [
            {
                'pdb': '1GID',
                'quality': {
                    'has': set([u'clashscore']),
                    'rsrz': 100,
                    'backbone': 100,
                    'clashscore': 36.87,
                }
            },
            {
                'pdb': '0S72',
                'quality': {
                    'has': set(),
                    'rsrz': 100,
                    'backbone': 100,
                    'clashscore': 500,
                }
            },
            {
                'pdb': '1FJG',
                'quality': {
                    'has': set(['rsrz', 'clashscore', 'backbone']),
                    'rsrz': 5.82,
                    'backbone': 7.0,
                    'clashscore': 29.48,
                }
            },
        ]


class FindingHardcodedTest(StageTest):
    loader_class = QualityMetrics

    def test_it_uses_exising_hardcoded(self):
        self.loader.hardcoded = set(['b', 'e', 'f'])
        assert self.loader.find_hardcoded([{'id': 'a'}, {'id': 'b'}]) == {'id': 'b'}

    def test_it_gives_none_for_no_hardcoded(self):
        self.loader.hardcoded = set()
        assert self.loader.find_hardcoded([{'id': 'a'}, {'id': 'b'}]) is None

    def test_it_gives_none_if_several_hardcoded(self):
        self.loader.hardcoded = set(['a', 'b', 'c'])
        assert self.loader.find_hardcoded([{'id': 'a'}, {'id': 'b'}]) is None


class UseHardcodedTest(StageTest):
    loader_class = QualityMetrics

    def test_it_places_hardcoded_first_if_exists(self):
        self.loader.hardcoded = set(['d', 'b', '0'])
        assert self.loader.use_hardcoded([{'id': 'a'}, {'id': 'b'}]) == [
            {'id': 'b'},
            {'id': 'a'}
        ]

    def test_it_does_not_alter_if_no_hardcoded(self):
        self.loader.hardcoded = set(['d', 'j', '0'])
        assert self.loader.use_hardcoded([{'id': 'a'}, {'id': 'b'}]) == [
            {'id': 'a'},
            {'id': 'b'}
        ]


class FilterByMethodTest(StageTest):
    loader_class = QualityMetrics

    def setUp(self):
        super(FilterByMethodTest, self).setUp()
        self.members = [
            {'id': 'a','method': 'X-RAY DIFFRACTION', 'resolution': 10, 'quality': {'has': set(['rsrz'])}},
            {'id': 'b','method': 'X-RAY DIFFRACTION', 'resolution': 4, 'quality': {'has': set(['rsrz', 'clashscore', 'backbone'])}},
            {'id': 'c','method': 'X-RAY DIFFRACTION', 'resolution': 2, 'quality': {'has': set([])}},
            {'id': 'c','method': 'X-RAY DIFFRACTION', 'resolution': 2, 'quality': {'has': set(['rsrz', 'clashscore', 'backbone'])}},
            {'id': 'd','method': 'X-RAY DIFFRACTION', 'resolution': 3, 'quality': {'has': set(['rsrz'])}},
            {'id': 'e','method': 'CRYO-EM', 'resolution': 2, 'quality': {'has': set(['rsrz'])}},
            {'id': 'f','method': 'CRYO-EM', 'resolution': 2, 'quality': {'has': set(['rsrz'])}},
        ]

    def filter(self, *args, **kwargs):
        return self.loader.filter_by_method(*args, **kwargs)

    def test_with_no_xray_uses_all(self):
        assert self.filter([self.members[6], self.members[5]]) == [
            self.members[6],
            self.members[5]
        ]

    def test_with_bad_xray_uses_all(self):
        assert self.filter([self.members[6], self.members[0], self.members[5]]) == [
            self.members[6],
            self.members[0],
            self.members[5]
        ]

    def test_with_xray_without_quality_uses_all(self):
        assert self.filter([self.members[2], self.members[5], self.members[6]]) == [
            self.members[2],
            self.members[5],
            self.members[6]
        ]

    def test_keeps_only_good_xray(self):
        assert self.filter([self.members[6], self.members[0], self.members[1], self.members[3], self.members[5]]) == [
            self.members[1],
            self.members[3],
        ]


class FilterByNtsTest(StageTest):
    loader_class = QualityMetrics

    def test_keeps_those_near_best(self):
        members = [
            {'id': 'a', 'length': 100},
            {'id': 'a', 'length': 90},
            {'id': 'a', 'length': 70},
            {'id': 'a', 'length': 80},
            {'id': 'a', 'length': 75},
        ]
        assert self.loader.filter_by_nts(members) == [
            members[0],
            members[1],
            members[3],
            members[4],
        ]


class FilterByResolutionTest(StageTest):
    loader_class = QualityMetrics

    def test_it_keeps_close_to_best(self):
        members = [
            {'id': 'a', 'resolution': 1},
            {'id': 'b', 'resolution': 1.1},
            {'id': 'c', 'resolution': 1.21},
            {'id': 'd', 'resolution': 1.2},
            {'id': 'e', 'resolution': 1.3},
        ]

        assert self.loader.filter_by_resolution(members) == [
            members[0],
            members[1],
            members[3],
        ]


class FinalOrderingTest(StageTest):
    loader_class = QualityMetrics

    def ordering(self, *args, **kwargs):
        members = self.loader.final_ordering(*args, **kwargs)
        return [m['id'] for m in members]

    def test_will_keep_given_order(self):
        assert self.ordering([{'id': 'a'}, {'id': 'b'}], []) == ['a', 'b']

    def test_will_reorder_remainder_by_quality(self):
        val = self.ordering(
            [{'id': 'a'}, {'id': 'b'}],
            [
                {
                    'id': 'e',
                    'resolution': 1,
                    'length': 10,
                    'bp': 10,
                    'quality': {
                        'rsrz': 100, 'backbone': 1000, 'clashscore': 1000,
                        'has': set(['rsrz', 'backbone', 'clashscore']),
                    },
                },
                {
                    'id': 'c',
                    'resolution': 1,
                    'length': 10,
                    'bp': 10,
                    'quality': {
                        'rsrz': 0, 'backbone': 0, 'clashscore': 0,
                        'has': set(['rsrz', 'backbone', 'clashscore']),
                    },
                },
            ])

        assert val == ['a', 'b', 'c', 'e']

    def test_does_not_list_ids_twice(self):
        val = self.ordering(
            [{'id': 'a'}, {'id': 'b'}],
            [
                {'id': 'a'},
                {'id': 'b'},
                {
                    'id': 'c',
                    'resolution': 1,
                    'length': 10,
                    'bp': 10,
                    'quality': {
                        'rsrz': 100, 'backbone': 1000, 'clashscore': 1000,
                        'has': set(['rsrz', 'backbone', 'clashscore']),
                    },
                },
                {
                    'resolution': 1,
                    'length': 10,
                    'bp': 10,
                    'id': 'e',
                    'quality': {
                        'rsrz': 0, 'backbone': 0, 'clashscore': 0,
                        'has': set(['rsrz', 'backbone', 'clashscore']),
                    },
                },
            ])

        assert val == ['a', 'b', 'e', 'c']


class RealDataTests(StageTest):
    loader_class = QualityMetrics

    def group(self, *args, **kwargs):
        grouper = Grouper(self.loader.config, self.loader.session)
        members = []
        for pdb, chains in args:
            ifes = grouper.ifes(pdb)
            members.extend(ife for ife in ifes if ife['name'] in chains)
        return members

    def rep(self, *args, **kwargs):
        group = self.group(*args, **kwargs)
        return [m['id'] for m in self.loader({'members': group})]

    def test_4MGM_4MGN(self):
        val = self.rep(('4MGM', ('A', 'B')), ('4MGN', ('B', 'D')),
                       id='NR_4.0_13428.1')
        assert val == ['4MGN|1|D', '4MGN|1|B', '4MGM|1|A', '4MGM|1|B']

    def test_1GID(self):
        val = self.rep(('1GID', ('A', 'B')), ('1X8W', ('A', 'B', 'C', 'D')),
                       ('1GRZ', ('A', 'B')), id='NR_4.0_86492.1')
        assert val == [
            '1X8W|1|B',
            '1X8W|1|A',
            '1X8W|1|D',
            '1X8W|1|C',
            '1GID|1|B',
            '1GID|1|A',
            '1GRZ|1|A',
            '1GRZ|1|B'
        ]
