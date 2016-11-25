import pytest

from test import StageTest

from pymotifs.core import InvalidState
from pymotifs.nr.builder import Known


class HandleTest(StageTest):
    loader_class = Known

    def test_can_get_all_known_handles(self):
        assert len(self.loader.handles()) == 88


class ClassesTest(StageTest):
    loader_class = Known

    def setUp(self):
        super(ClassesTest, self).setUp()
        self.data = {}
        for klass in self.loader.classes('1.0', '4.0'):
            self.data[klass['name']['full']] = klass

    def test_can_get_all_known_classes(self):
        assert len(self.data) == 31

    def test_can_compute_correct_data(self):
        assert self.data['NR_4.0_13436.1'] == {
            'members': [
                {'id': '4V88|1|A7'},
                {'id': '4V88|1|A3'},
                {'id': '4V7R|1|D2'},
                {'id': '4V7R|1|B2'}
            ],
            'representative': {'id': '4V88|1|A7', 'length': 121, 'bp': 53},
            'name': {
                'class_id': 97,
                'full': 'NR_4.0_13436.1',
                'handle': '13436',
                'cutoff': '4.0',
                'version': 1
            },
            'release': '1.0'
        }


class MappingTest(StageTest):
    loader_class = Known

    def test_can_generate_mapping_for_names(self):
        assert self.loader.mapping('1.0', ['NR_4.0_13436.1']) == {
            'NR_4.0_13436.1': 97
        }

    def test_complains_if_no_release(self):
        with pytest.raises(InvalidState):
            self.loader.mapping(None, ['a'])

    def test_complains_if_no_names(self):
        with pytest.raises(InvalidState):
            self.loader.mapping('1.0', [])

    def test_complains_if_cannot_map_all_names(self):
        with pytest.raises(InvalidState):
            self.loader.mapping('1.0', ['a', 'NR_4.0_13436.1'])

    def test_complains_if_finds_no_matches(self):
        with pytest.raises(InvalidState):
            self.loader.mapping('1.0', ['a'])
