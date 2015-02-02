from nose import SkipTest

from test import QueryUtilTest

from pymotifs.utils import structures as st


class RepresentativeBpTest(QueryUtilTest):
    query_class = st.BasePairQueries

    def test_can_get_no_interactions(self):
        raise SkipTest("No example data yet")

    def test_gets_count_in_a_simple_structure(self):
        val = self.db_obj.representative('1EIY', 'C', count=True)
        ans = 27
        self.assertEquals(ans, val)

    def test_gets_count_for_one_chain(self):
        val = self.db_obj.representative('1J5E', 'A', count=True)
        ans = 689
        self.assertEquals(ans, val)

    def test_gets_count_for_several_chains(self):
        val = self.db_obj.representative('3U5H', ['5', '8'], count=True)
        ans = 1321
        self.assertEquals(ans, val)

    def test_it_can_get_the_long_range_count(self):
        val = self.db_obj.representative('1EIY', 'C', count=True,
                                         range_cutoff=st.LONG_RANGE_CUTOFF)
        ans = 3
        self.assertEquals(ans, val)


class CrossChainTest(QueryUtilTest):
    query_class = st.BasePairQueries

    def test_gets_zero_count(self):
        val = self.db_obj.cross_chain('1J5E', 'A', count=True)
        ans = 0
        self.assertEquals(ans, val)

    def test_counts_from_a_chosen_loop(self):
        val = self.db_obj.cross_chain('1ET4', 'A', count=True)
        ans = 3
        self.assertEquals(ans, val)

    def test_counts_from_several_chosen_loops(self):
        val = self.db_obj.cross_chain('1ET4', ['A', 'B'], count=True)
        ans = 4
        self.assertEquals(ans, val)

    def test_counts_between_specific_chains(self):
        val = self.db_obj.cross_chain('1ET4', 'A', other_chain='B',
                                      count=True)
        ans = 3
        self.assertEquals(ans, val)

    def test_can_get_a_zero_count_to_other_chains(self):
        val = self.db_obj.cross_chain('1ET4', 'A', other_chain='D',
                                      count=True)
        ans = 0
        self.assertEquals(ans, val)

    def test_counts_to_a_list_of_other_chains(self):
        val = self.db_obj.cross_chain('1ET4', 'A', other_chain=['B', 'E'],
                                      count=True)
        ans = 3
        self.assertEquals(ans, val)


class BetweenChainTest(QueryUtilTest):
    query_class = st.BasePairQueries

    def test_it_can_get_a_zero_count(self):
        val = self.db_obj.between('1J5E', ['A', 'C'], count=True)
        ans = 0
        self.assertEquals(ans, val)

    def test_counts_from_several_chosen_loops(self):
        val = self.db_obj.cross_chain('1ET4', ['A', 'B'], count=True)
        ans = 4
        self.assertEquals(ans, val)
