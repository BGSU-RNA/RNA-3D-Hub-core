import pytest

from test import StageTest

from pymotifs.nr.parents import Loader


class ParentInfoTest(StageTest):
    loader_class = Loader

    def test_build_valid_data(self):
        val = self.loader.parents('1.0', [
            {'name': {'class_id': 1}, 'parents': [{'name': {'class_id': -2}},
                                                  {'name': {'class_id': -3}}]},
            {'name': {'class_id': 3}, 'parents': []},
        ])
        assert val == [
            {'nr_class_id': 1, 'nr_release_id': '1.0', 'nr_class_parent_id': -2},
            {'nr_class_id': 1, 'nr_release_id': '1.0', 'nr_class_parent_id': -3},
        ]

    def test_knows_if_no_parents(self):
        data = {'parent_counts': [
            {'cutoff': '1.5', 'classes': {'updated': 0, 'removed': 61, 'unchanged': 0, 'added': 65}, 'pdbs': {'removed': 87, 'added': 1}, 'ifes': {'unchanged': 0, 'removed': 173, 'added': 197}},
            {'cutoff': '2.0', 'classes': {'updated': 0, 'removed': 224, 'unchanged': 0, 'added': 242}, 'pdbs': {'removed': 226, 'added': 6}, 'ifes': {'unchanged': 0, 'removed': 552, 'added': 587}},
            {'cutoff': '2.5', 'classes': {'updated': 0, 'removed': 447, 'unchanged': 0, 'added': 507}, 'pdbs': {'removed': 786, 'added': 8}, 'ifes': {'unchanged': 0, 'removed': 1499, 'added': 1234}},
            {'cutoff': '20.0', 'classes': {'updated': 0, 'removed': 945, 'unchanged': 0, 'added': 1353}, 'pdbs': {'removed': 760, 'added': 272}, 'ifes': {'unchanged': 0, 'removed': 2588, 'added': 4931}},
            {'cutoff': '3.0', 'classes': {'updated': 0, 'removed': 700, 'unchanged': 0, 'added': 803}, 'pdbs': {'removed': 964, 'added': 50}, 'ifes': {'unchanged': 0, 'removed': 2249, 'added': 2609}},
            {'cutoff': '3.5', 'classes': {'updated': 0, 'removed': 838, 'unchanged': 0, 'added': 1007}, 'pdbs': {'removed': 814, 'added': 125}, 'ifes': {'unchanged': 0, 'removed': 2430, 'added': 3754}},
            {'cutoff': '4.0', 'classes': {'updated': 0, 'removed': 876, 'unchanged': 0, 'added': 1075}, 'pdbs': {'removed': 766, 'added': 160}, 'ifes': {'unchanged': 0, 'removed': 2488, 'added': 4143}},
            {'cutoff': 'all', 'classes': {'updated': 0, 'removed': 1419, 'unchanged': 0, 'added': 1788}, 'pdbs': {'removed': 737, 'added': 272}, 'ifes': {'unchanged': 0, 'removed': 3145, 'added': 5549}}
        ]}
        assert self.loader.no_parents(data) is True

    def test_knows_if_has_parents(self):
        updated = {'parent_counts': [
            {'cutoff': '1.5', 'classes': {'updated': 1, 'removed': 61, 'unchanged': 0, 'added': 65}, 'pdbs': {'removed': 87, 'added': 1}, 'ifes': {'unchanged': 0, 'removed': 173, 'added': 197}},
            {'cutoff': '2.0', 'classes': {'updated': 1, 'removed': 224, 'unchanged': 0, 'added': 242}, 'pdbs': {'removed': 226, 'added': 6}, 'ifes': {'unchanged': 0, 'removed': 552, 'added': 587}},
            {'cutoff': '2.5', 'classes': {'updated': 1, 'removed': 447, 'unchanged': 0, 'added': 507}, 'pdbs': {'removed': 786, 'added': 8}, 'ifes': {'unchanged': 0, 'removed': 1499, 'added': 1234}},
            {'cutoff': '20.0', 'classes': {'updated': 1, 'removed': 945, 'unchanged': 0, 'added': 1353}, 'pdbs': {'removed': 760, 'added': 272}, 'ifes': {'unchanged': 0, 'removed': 2588, 'added': 4931}},
            {'cutoff': '3.0', 'classes': {'updated': 1, 'removed': 700, 'unchanged': 0, 'added': 803}, 'pdbs': {'removed': 964, 'added': 50}, 'ifes': {'unchanged': 0, 'removed': 2249, 'added': 2609}},
            {'cutoff': '3.5', 'classes': {'updated': 1, 'removed': 838, 'unchanged': 0, 'added': 1007}, 'pdbs': {'removed': 814, 'added': 125}, 'ifes': {'unchanged': 0, 'removed': 2430, 'added': 3754}},
            {'cutoff': '4.0', 'classes': {'updated': 1, 'removed': 876, 'unchanged': 0, 'added': 1075}, 'pdbs': {'removed': 766, 'added': 160}, 'ifes': {'unchanged': 0, 'removed': 2488, 'added': 4143}},
            {'cutoff': 'all', 'classes': {'updated': 1, 'removed': 1419, 'unchanged': 0, 'added': 1788}, 'pdbs': {'removed': 737, 'added': 272}, 'ifes': {'unchanged': 0, 'removed': 3145, 'added': 5549}}
        ]}
        unchanged = {'parent_counts': [
            {'cutoff': '1.5', 'classes': {'updated': 0, 'removed': 61, 'unchanged': 1, 'added': 65}, 'pdbs': {'removed': 87, 'added': 1}, 'ifes': {'unchanged': 0, 'removed': 173, 'added': 197}},
            {'cutoff': '2.0', 'classes': {'updated': 0, 'removed': 224, 'unchanged': 1, 'added': 242}, 'pdbs': {'removed': 226, 'added': 6}, 'ifes': {'unchanged': 0, 'removed': 552, 'added': 587}},
            {'cutoff': '2.5', 'classes': {'updated': 0, 'removed': 447, 'unchanged': 1, 'added': 507}, 'pdbs': {'removed': 786, 'added': 8}, 'ifes': {'unchanged': 0, 'removed': 1499, 'added': 1234}},
            {'cutoff': '20.0', 'classes': {'updated': 0, 'removed': 945, 'unchanged': 1, 'added': 1353}, 'pdbs': {'removed': 760, 'added': 272}, 'ifes': {'unchanged': 0, 'removed': 2588, 'added': 4931}},
            {'cutoff': '3.0', 'classes': {'updated': 0, 'removed': 700, 'unchanged': 1, 'added': 803}, 'pdbs': {'removed': 964, 'added': 50}, 'ifes': {'unchanged': 0, 'removed': 2249, 'added': 2609}},
            {'cutoff': '3.5', 'classes': {'updated': 0, 'removed': 838, 'unchanged': 1, 'added': 1007}, 'pdbs': {'removed': 814, 'added': 125}, 'ifes': {'unchanged': 0, 'removed': 2430, 'added': 3754}},
            {'cutoff': '4.0', 'classes': {'updated': 0, 'removed': 876, 'unchanged': 1, 'added': 1075}, 'pdbs': {'removed': 766, 'added': 160}, 'ifes': {'unchanged': 0, 'removed': 2488, 'added': 4143}},
            {'cutoff': 'all', 'classes': {'updated': 0, 'removed': 1419, 'unchanged': 1, 'added': 1788}, 'pdbs': {'removed': 737, 'added': 272}, 'ifes': {'unchanged': 0, 'removed': 3145, 'added': 5549}}
        ]}
        both = {'parent_counts': [
            {'cutoff': '1.5', 'classes': {'updated': 1, 'removed': 61, 'unchanged': 1, 'added': 65}, 'pdbs': {'removed': 87, 'added': 1}, 'ifes': {'unchanged': 0, 'removed': 173, 'added': 197}},
            {'cutoff': '2.0', 'classes': {'updated': 1, 'removed': 224, 'unchanged': 1, 'added': 242}, 'pdbs': {'removed': 226, 'added': 6}, 'ifes': {'unchanged': 0, 'removed': 552, 'added': 587}},
            {'cutoff': '2.5', 'classes': {'updated': 1, 'removed': 447, 'unchanged': 1, 'added': 507}, 'pdbs': {'removed': 786, 'added': 8}, 'ifes': {'unchanged': 0, 'removed': 1499, 'added': 1234}},
            {'cutoff': '20.0', 'classes': {'updated': 2, 'removed': 945, 'unchanged': 1, 'added': 1353}, 'pdbs': {'removed': 760, 'added': 272}, 'ifes': {'unchanged': 0, 'removed': 2588, 'added': 4931}},
            {'cutoff': '3.0', 'classes': {'updated': 1, 'removed': 700, 'unchanged': 1, 'added': 803}, 'pdbs': {'removed': 964, 'added': 50}, 'ifes': {'unchanged': 0, 'removed': 2249, 'added': 2609}},
            {'cutoff': '3.5', 'classes': {'updated': 1, 'removed': 838, 'unchanged': 1, 'added': 1007}, 'pdbs': {'removed': 814, 'added': 125}, 'ifes': {'unchanged': 0, 'removed': 2430, 'added': 3754}},
            {'cutoff': '4.0', 'classes': {'updated': 1, 'removed': 876, 'unchanged': 1, 'added': 1075}, 'pdbs': {'removed': 766, 'added': 160}, 'ifes': {'unchanged': 0, 'removed': 2488, 'added': 4143}},
            {'cutoff': 'all', 'classes': {'updated': 1, 'removed': 1419, 'unchanged': 1, 'added': 1788}, 'pdbs': {'removed': 737, 'added': 272}, 'ifes': {'unchanged': 0, 'removed': 3145, 'added': 5549}}
        ]}
        assert self.loader.no_parents(updated) is False
        assert self.loader.no_parents(unchanged) is False
        assert self.loader.no_parents(both) is False

    @pytest.mark.skip()
    def test_it_gets_all_parent_names(self):
        val, _ = self.loader.parents(self.grouping)
        ans = sorted(['NR_1.5_01181.1', 'NR_1.5_08345.1', 'NR_1.5_23793.1',
                      'NR_1.5_26877.1', 'NR_1.5_31163.1'])
        self.assertEquals(ans, sorted(val))

    @pytest.mark.skip()
    def test_gives_no_parents_if_no_previous_release(self):
        val = self.loader.parents(self.grouping)
        ans = ([], None)
        self.assertEquals(ans, val)
