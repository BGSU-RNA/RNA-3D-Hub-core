from test import StageTest

from pymotifs.utils.naming import ChangeCounter

GROUPS = [
    {'members': [{'id': 1}, {'id': 2}, {'id': -2}], 'name': {'handle': 1, 'version': 2}},
    {'members': [{'id': 3}, {'id': 4}, {'id': -1}], 'name': {'handle': 2, 'version': 1}},
    {'members': [{'id': 5}, {'id': 6}], 'name': {'handle': 4, 'version': 1}},
    {'members': [{'id': 7}, {'id': 8}], 'name': {'handle': 5, 'version': 1}},
    {'members': [{'id': 9}, {'id': 10}], 'name': {'handle': 6, 'version': 1}},
    {'members': [{'id': -9}, {'id': -10}], 'name': {'handle': 'a', 'version': 1}},
    {'members': [{'id': -7}, {'id': -8}], 'name': {'handle': 'b', 'version': 1}},
    {'members': [{'id': -5}, {'id': -6}], 'name': {'handle': 'c', 'version': 1}},
]

PARENTS = [
    {'members': [{'id': 1}, {'id': 'a'}], 'name': {'handle': 1, 'version': 1}},
    {'members': [{'id': 3}, {'id': 'b'}], 'name': {'handle': 3, 'version': 1}},
    {'members': [{'id': 5}, {'id': 6}], 'name': {'handle': 4, 'version': 1}},
    {'members': [{'id': 7}, {'id': 8}], 'name': {'handle': 5, 'version': 1}},
    {'members': [{'id': 9}, {'id': 10}], 'name': {'handle': 6, 'version': 1}},
]


class ChangesTest(StageTest):
    loader_class = ChangeCounter

    def test_can_get_correct_number_of_group_changes(self):
        assert self.loader.group_changes(GROUPS, PARENTS) == {
            'added': [GROUPS[1], GROUPS[5], GROUPS[6], GROUPS[7]],
            'removed': [PARENTS[1]],
            'updated': [GROUPS[0]],
            'unchanged': GROUPS[2:5],
        }

    def test_can_get_correct_number_of_transformed_counts(self):
        def fn(group):
            return [m['id'] for m in group['members']]

        val = self.loader.transformed_changes(GROUPS, PARENTS, fn)
        assert val == {
            'added': set([-2, -1, -10, -9, -8, -7, -6, -5, 2, 4]),
            'removed': set(['a', 'b']),
            'unchanged': set([1, 3, 5, 6, 7, 8, 9, 10]),
        }


class CountsTests(StageTest):
    loader_class = ChangeCounter

    def test_can_get_correct_counts(self):
        def fn(group):
            val = []
            for member in group['members']:
                if member['id'] > 0:
                    val.append(-1 * member['id'])
                else:
                    val.append(member['id'])
            return val

        assert self.loader(GROUPS, PARENTS, example=fn) == {
            'groups': {
                'added': 4,
                'removed': 1,
                'updated': 1,
                'unchanged': 3,
            },
            'members': {
                'added': 10,
                'removed': 2,
                'unchanged': 8,
            },
            'example': {
                'added': 2,
                'removed': 1,
                'unchanged': 8,
            }
        }
