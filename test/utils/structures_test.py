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

    def test_it_can_limit_by_family(self):
        val = self.db_obj.representative('1J5E', 'A', count=True, family='tHH')
        ans = 4
        self.assertEquals(ans, val)

    def test_it_can_get_the_long_range_count(self):
        val = self.db_obj.representative('1EIY', 'C', count=True,
                                         range_cutoff=st.LONG_RANGE_CUTOFF)
        ans = 3
        self.assertEquals(ans, val)

    def test_it_counts_only_in_one_model(self):
        raise SkipTest()

    def test_it_counts_only_in_one_symmetry_operator(self):
        raise SkipTest()


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

    def test_it_can_limit_by_family(self):
        val = self.db_obj.cross_chain('1ET4', 'A', other_chain=['B', 'E'],
                                      count=True, family='tWS')
        ans = 1
        self.assertEquals(ans, val)

    def test_it_can_find_only_in_one_model(self):
        raise SkipTest()

    def test_it_can_find_in_only_one_symmetry_operator(self):
        raise SkipTest()


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

    def test_can_limit_by_family(self):
        val = self.db_obj.cross_chain('1ET4', ['A', 'B'], count=True,
                                      family='tWS')
        ans = 2
        self.assertEquals(ans, val)

    def test_it_searches_only_in_one_model(self):
        raise SkipTest()

    def test_it_searches_only_in_one_symmetry_operator(self):
        raise SkipTest()


class SourceTest(QueryUtilTest):
    query_class = st.Structure

    def test_can_find_a_simple_taxon_id(self):
        self.assertEquals([562], self.db_obj.source('2AW7', 'A'))

    def test_knows_when_there_is_no_known_taxon_id(self):
        self.assertEquals([], self.db_obj.source('1ET4', 'A'))

    def test_can_find_several_species_ids(self):
        self.assertEquals([11103, 32630], self.db_obj.source('3T4B', 'A'))

    def test_fails_if_it_cannot_find_all_taxon_ids(self):
        raise SkipTest()

    def test_can_find_by_chain_id(self):
        self.assertEquals([562], self.db_obj.source('2AW7', 4977))

    def test_can_simplify_to_just_syntheic(self):
        self.assertEquals(32630, self.db_obj.source('3T4B', 'A',
                                                    simplify=True))

    def test_simplifies_to_just_first_id(self):
        self.assertEquals(562, self.db_obj.source('2AW7', 'A', simplify=True))


class RnaChainsTest(QueryUtilTest):
    query_class = st.Structure

    def test_can_get_all_rna_chains_in_a_structure(self):
        self.assertEquals(['A'], self.db_obj.rna_chains('3T4B'))

    def test_will_filter_out_non_rna_chains(self):
        self.assertEquals(['A'], self.db_obj.rna_chains('2AW7'))

    def test_can_return_ids_and_names(self):
        val = self.db_obj.rna_chains('2AW7', return_id=True)
        self.assertEquals([('A', 4977)], val)
