import pytest

from test import QueryUtilTest

from pymotifs.utils import structures as st


class RepresentativeBpTest(QueryUtilTest):
    query_class = st.BasePairQueries

    @pytest.mark.skip("No example data yet")
    def test_can_get_no_interactions(self):
        pass

    def test_gets_count_in_a_simple_structure(self):
        val = self.db_obj.representative('1EIY', 'C', count=True)
        ans = 27
        self.assertEquals(ans, val)

    def test_gets_count_for_one_chain(self):
        val = self.db_obj.representative('1J5E', 'A', count=True)
        ans = 691
        self.assertEquals(ans, val)

    def test_gets_count_for_several_chains(self):
        val = self.db_obj.representative('1GID', ['A', 'B'], count=True)
        ans = 149
        self.assertEquals(ans, val)

    def test_it_can_limit_by_family(self):
        val = self.db_obj.representative('1J5E', 'A', count=True, family='tHH')
        ans = 4
        self.assertEquals(ans, val)

    def test_it_can_get_the_long_range_count(self):
        val = self.db_obj.representative('1EIY', 'C', count=True,
                                         range_cutoff=st.LONG_RANGE_CUTOFF)
        ans = 3
        self.assertEquals(ans, val)

    def test_can_load_cWW_given_symmetry_ops(self):
        val = self.db_obj.representative('4PMI', 'A', count=True, family='cWW',
                                         sym_op='6_445')
        self.assertEquals(val, 15)

    def test_can_load_all_given_symmetry_ops(self):
        val = self.db_obj.representative('4PMI', 'A', count=True,
                                         sym_op='6_445')
        self.assertEquals(val, 17)

    @pytest.mark.skip("No example data yet")
    def test_it_counts_only_in_one_model(self):
        pass

    def test_it_can_count_whole_structure(self):
        val = self.db_obj.representative('4PMI', None, count=True)
        self.assertEquals(val, 17)


class CrossChainTest(QueryUtilTest):
    query_class = st.BasePairQueries

    def test_gets_zero_count(self):
        val = self.db_obj.cross_chain('1J5E', 'A', count=True)
        ans = 0
        self.assertEquals(ans, val)

    def test_counts_from_a_chosen_chain(self):
        val = self.db_obj.cross_chain('1ET4', 'A', count=True)
        ans = 3
        self.assertEquals(ans, val)

    def test_counts_from_several_chosen_chains(self):
        val = self.db_obj.cross_chain('1ET4', ['A', 'B'], count=True)
        ans = 6
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

    def test_it_can_limit_by_family(self):
        val = self.db_obj.cross_chain('1ET4', 'A', other_chain=['B', 'E'],
                                      count=True, family='tWS')
        ans = 1
        self.assertEquals(ans, val)

    @pytest.mark.skip("No example data yet")
    def test_it_can_find_only_in_one_model(self):
        pass

    @pytest.mark.skip("No example data yet")
    def test_it_can_find_in_only_one_symmetry_operator(self):
        pass


class BetweenChainTest(QueryUtilTest):
    query_class = st.BasePairQueries

    def test_it_can_get_a_zero_count(self):
        val = self.db_obj.between('1J5E', ['A', 'C'], count=True)
        ans = 0
        self.assertEquals(ans, val)

    def test_counts_from_several_chosen_chains(self):
        val = self.db_obj.cross_chain('1ET4', ['A', 'B'], count=True)
        ans = 6
        self.assertEquals(ans, val)

    def test_can_limit_by_family(self):
        val = self.db_obj.cross_chain('1ET4', ['A', 'B'], count=True,
                                      family='tWS')
        ans = 2
        self.assertEquals(ans, val)

    @pytest.mark.skip("No example data yet")
    def test_it_searches_only_in_one_model(self):
        pass

    @pytest.mark.skip("No example data yet")
    def test_it_searches_only_in_one_symmetry_operator(self):
        pass
