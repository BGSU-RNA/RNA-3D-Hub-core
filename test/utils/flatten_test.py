from pymotifs.utils import flatten

from unittest import TestCase


class FlattenTests(TestCase):
    def test_it_can_flatten_a_list(self):
        val = list(flatten([1, [[2], []], 3, (4,)]))
        assert val == [1, 2, 3, 4]
