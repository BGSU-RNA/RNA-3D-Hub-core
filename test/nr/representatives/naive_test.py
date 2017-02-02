from test import StageTest

from pymotifs.nr.representatives import Naive


class NaiveBestTest(StageTest):
    loader_class = Naive

    def best(self, *chains):
        ordered = self.loader({'members': chains, 'parent': []})
        return ordered[0]

    def test_it_can_deal_with_0_bp_or_length_chains(self):
        val = self.best({'bp': 13, 'length': 10, 'id': 'c'},
                        {'bp': 0, 'length': 10, 'id': 'a'},
                        {'bp': 50, 'length': 0, 'id': 'b'})
        assert val == {'bp': 13, 'length': 10, 'id': 'c'}

    def test_gets_by_bp_per_nt(self):
        val = self.best({'bp': 13, 'length': 10, 'id': 'c'},
                        {'bp': 10, 'length': 10, 'id': 'a'},
                        {'bp': 50, 'length': 10, 'id': 'b'})
        ans = {'bp': 50, 'length': 10, 'id': 'b'}
        self.assertEquals(ans, val)

    def test_it_tiebreaks_on_pdb(self):
        val = self.best({'bp': 10, 'length': 10, 'id': 'c'},
                        {'bp': 10, 'length': 10, 'id': 'a'},
                        {'bp': 10, 'length': 10, 'id': 'b'})
        ans = {'bp': 10, 'length': 10, 'id': 'c'}
        self.assertEquals(ans, val)

    def test_if_all_have_0_bp_uses_resolution(self):
        val = self.best({'bp': 0, 'length': 10, 'id': 'c', 'resolution': 4.0},
                        {'bp': 0, 'length': 10, 'id': 'a', 'resolution': 3.0},
                        {'bp': 0, 'length': 10, 'id': 'b', 'resolution': 2.0})
        ans = {'bp': 0, 'length': 10, 'id': 'b', 'resolution': 2.0}
        self.assertEquals(ans, val)

    def test_it_ignores_any_set_parent(self):
        val = self.loader({
            'parent': [{'representative': {'bp': 100, 'length': 10}}],
            'members': [{'bp': 13, 'length': 10, 'id': 'c'},
                        {'bp': 0, 'length': 10, 'id': 'a'},
                        {'bp': 50, 'length': 0, 'id': 'b'}],
        })
        assert val == [
            {'bp': 13, 'length': 10, 'id': 'c'},
            {'bp': 50, 'length': 0, 'id': 'b'},
            {'bp': 0, 'length': 10, 'id': 'a'},
        ]

    def test_it_orders_by_bp_nt(self):
        ordered = self.loader({
            'members': [{'bp': 13, 'length': 10, 'id': 'c'},
                        {'bp': 0, 'length': 10, 'id': 'a'},
                        {'bp': 50, 'length': 0, 'id': 'b'}]
        })
        assert [m['id'] for m in ordered] == ['c', 'b', 'a']
