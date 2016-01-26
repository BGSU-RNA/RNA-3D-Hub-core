from test import StageTest
from nose import SkipTest

from pymotifs.nr.builder import RepresentativeFinder


class SortingChainsTest(StageTest):
    loader_class = RepresentativeFinder

    def key(self, chain):
        return self.loader.sorting_key(chain)

    def test_builds_key_with_nt_per_bp_and_id(self):
        chain1 = {'bps': 10, 'length': 20, 'id': 'a'}
        chain2 = {'bps': 20, 'length': 20, 'id': 'b'}
        self.assertEquals(self.key(chain1), (0.5, 'a'))
        self.assertEquals(self.key(chain2), (1, 'b'))

    def test_builds_key_with_0_if_one_is_0(self):
        chain1 = {'bps': 0, 'length': 10, 'id': 'c'}
        chain2 = {'bps': 10, 'length': 0, 'id': 'd'}
        chain3 = {'bps': 0, 'length': 0, 'id': 'e'}
        self.assertEquals(self.key(chain1), (0, 'c'))
        self.assertEquals(self.key(chain2), (0, 'd'))
        self.assertEquals(self.key(chain3), (0, 'e'))


class NaiveBestTest(StageTest):
    loader_class = RepresentativeFinder

    def best(self, *chains):
        return self.loader.naive_best(chains)

    def test_it_can_deal_with_0_bp_or_length_chains(self):
        val = self.best({'bps': 13, 'length': 10, 'id': 'c'},
                        {'bps': 0, 'length': 10, 'id': 'a'},
                        {'bps': 50, 'length': 0, 'id': 'b'})
        ans = {'bps': 13, 'length': 10, 'id': 'c'}
        self.assertEquals(ans, val)

    def test_gets_by_bps_per_nt(self):
        val = self.best({'bps': 13, 'length': 10, 'id': 'c'},
                        {'bps': 10, 'length': 10, 'id': 'a'},
                        {'bps': 50, 'length': 10, 'id': 'b'})
        ans = {'bps': 50, 'length': 10, 'id': 'b'}
        self.assertEquals(ans, val)

    def test_it_tiebreaks_on_pdb(self):
        val = self.best({'bps': 10, 'length': 10, 'id': 'c'},
                        {'bps': 10, 'length': 10, 'id': 'a'},
                        {'bps': 10, 'length': 10, 'id': 'b'})
        ans = {'bps': 10, 'length': 10, 'id': 'c'}
        self.assertEquals(ans, val)


class CandidatesTest(StageTest):
    loader_class = RepresentativeFinder

    def candidates(self, *chains):
        default = {'bps': 10, 'length': 20, 'id': 'a'}
        return self.loader.candidates(default, chains)

    def test_it_can_get_all_longer_ifes(self):
        val = self.candidates({'bps': 10, 'length': 10, 'id': 'b'},
                              {'bps': 9, 'length': 100, 'id': 'e'},
                              {'bps': 10, 'length': 30, 'id': 'c'},
                              {'bps': 10, 'length': 20, 'id': 'd'})
        ans = [{'bps': 10, 'length': 30, 'id': 'c'}]
        self.assertEquals(ans, val)

    def test_it_can_get_all_with_more_bps(self):
        val = self.candidates({'bps': 10, 'length': 20, 'id': 'b'},
                              {'bps': 100, 'length': 1, 'id': 'e'},
                              {'bps': 13, 'length': 20, 'id': 'c'},
                              {'bps': 0, 'length': 20, 'id': 'd'})
        ans = [{'bps': 13, 'length': 20, 'id': 'c'}]
        self.assertEquals(ans, val)

    def test_it_can_sort_results_correctly(self):
        val = self.candidates({'bps': 11, 'length': 20, 'id': 'b'},
                              {'bps': 0, 'length': 20, 'id': 'd'},
                              {'bps': 13, 'length': 20, 'id': 'c'})
        ans = [{'bps': 11, 'length': 20, 'id': 'b'},
               {'bps': 13, 'length': 20, 'id': 'c'}]
        self.assertEquals(ans, val)


class IncreaseTest(StageTest):
    loader_class = RepresentativeFinder

    def increase(self, chain1, chain2, key='a'):
        return self.loader.increase(chain1, chain2, key)

    def test_it_computes_percent_increase(self):
        val = self.increase({'a': 10}, {'a': 1})
        self.assertEquals(val, 900.0)

    def test_it_can_handle_stays_same(self):
        val = self.increase({'a': 10}, {'a': 10})
        self.assertEquals(val, 0.0)

    def test_it_can_handle_decrease(self):
        val = self.increase({'a': 5}, {'a': 10})
        self.assertEquals(val, -50.0)

    def test_it_gives_100_for_given_0(self):
        self.assertEquals(self.increase({'a': 0}, {'a': 10}), -100)
        self.assertEquals(self.increase({'a': 10}, {'a': 0}), 100)

    def test_given_both_0_it_gives_0(self):
        self.assertEquals(self.increase({'a': 0}, {'a': 0}), 0)


class BestAboveCutoffsTest(StageTest):
    loader_class = RepresentativeFinder

    def test_it_finds_chain_with_more_bps_and_nts(self):
        raise SkipTest()

    def test_it_will_increase_if_length_increases(self):
        raise SkipTest()

    def test_it_will_increase_if_bps_increase(self):
        raise SkipTest()


class PickingRepresentativeTest(StageTest):
    loader_class = RepresentativeFinder
