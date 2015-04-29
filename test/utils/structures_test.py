from nose import SkipTest

from test import QueryUtilTest

from pymotifs.utils import structures as st


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

    def test_does_not_load_mislabeled_chain(self):
        val = self.db_obj.rna_chains('3CPW', strict=True)
        self.assertTrue('A' not in val)


class LoopLoadingTest(QueryUtilTest):
    query_class = st.Structure

    def test_can_load_data(self):
        raise SkipTest
        val = self.db_obj.loops('IL_1FG0_007')
        ans = {
            'id': 'IL_1FG0_007',
            'nts': [],
            'endpoints': [('1FG0_AU_1_A_2331_C_', '1FG0_AU_1_A_2333_G_'),
                          ('1FG0_AU_1_A_2351_C_', '1FG0_AU_1_A_2355_G_')]
        }
        self.assertEquals(ans, val)
