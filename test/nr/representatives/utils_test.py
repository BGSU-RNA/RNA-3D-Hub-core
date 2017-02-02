from test import StageTest

from pymotifs.nr.representatives import bp_per_nt


class SortingChainsTest(StageTest):

    def test_builds_key_with_nt_per_bp_and_id(self):
        chain1 = {'bp': 10, 'length': 20, 'id': 'a'}
        chain2 = {'bp': 20, 'length': 20, 'id': 'b'}
        self.assertEquals(bp_per_nt(chain1), (0.5, None, 'a'))
        self.assertEquals(bp_per_nt(chain2), (1, None, 'b'))

    def test_builds_key_with_0_if_one_is_0(self):
        chain1 = {'bp': 0, 'length': 10, 'id': 'c'}
        chain2 = {'bp': 10, 'length': 0, 'id': 'd'}
        chain3 = {'bp': 0, 'length': 0, 'id': 'e'}
        self.assertEquals(bp_per_nt(chain1), (0, None, 'c'))
        self.assertEquals(bp_per_nt(chain2), (0, None, 'd'))
        self.assertEquals(bp_per_nt(chain3), (0, None, 'e'))

    def test_uses_inverse_resolution(self):
        chain1 = {'bp': 0, 'length': 10, 'id': 'c', 'resolution': 10}
        self.assertEquals(bp_per_nt(chain1), (0, -10, 'c'))
