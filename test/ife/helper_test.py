from unittest import TestCase

from test import StageTest

from pymotifs.ife.helpers import IfeChain
from pymotifs.ife.helpers import IfeGroup
from pymotifs.ife.helpers import IfeLoader


class IfeChainTest(TestCase):
    def test_knows_if_chain_is_structured(self):
        val = IfeChain(internal=5)
        self.assertTrue(val.is_structured)

    def test_knows_if_chain_is_not_structured(self):
        val = IfeChain(internal=4)
        self.assertFalse(val.is_structured)

    def test_knows_if_two_chains_are_equal(self):
        val1 = IfeChain(internal=4)
        val2 = IfeChain(internal=4)
        self.assertEquals(val1, val2)

    def test_knows_if_two_chains_are_not_equal(self):
        val1 = IfeChain(internal=4)
        val2 = IfeChain(chain='G', internal=4)
        self.assertNotEquals(val1, val2)

    def test_sorts_by_chain(self):
        val1 = IfeChain(chain='A')
        val2 = IfeChain(chain='B')
        self.assertTrue(val1 > val2)
        self.assertTrue(val2 < val1)

    def test_sorts_by_length(self):
        val1 = IfeChain(pdb='01GG', chain='A', internal=0, length=10,
                        full_length=10, db_id=1)
        val2 = IfeChain(pdb='01GG', chain='A', internal=0, length=12,
                        full_length=10, db_id=1)
        self.assertTrue(val1 < val2)
        self.assertTrue(val2 > val1)

    def test_sorts_by_full_length(self):
        val1 = IfeChain(pdb='01GG', chain='A', internal=0, length=10,
                        full_length=10, db_id=1)
        val2 = IfeChain(pdb='01GG', chain='A', internal=0, length=10,
                        full_length=12, db_id=1)
        self.assertTrue(val1 < val2)
        self.assertTrue(val2 > val1)

    def test_it_is_bigger_than_none(self):
        val = IfeChain(pdb='01GG', chain='A', internal=0, length=10,
                       full_length=10, db_id=1)
        self.assertTrue(val > None)
        self.assertTrue(None < val)

    def test_compares_using_length(self):
        ife1 = IfeChain(pdb='0111', chain='A', length=10, internal=50)
        ife2 = IfeChain(pdb='0111', chain='C', length=50, internal=50)
        self.assertTrue(ife1 < ife2)
        self.assertTrue(ife1 <= ife2)
        self.assertFalse(ife1 == ife2)
        self.assertTrue(ife2 > ife1)

    def test_it_sorts_by_chain_and_internal(self):
        ife1 = IfeChain(pdb='0111', chain='A', internal=5)
        ife2 = IfeChain(pdb='0111', chain='D', internal=6)
        self.assertTrue(ife1 < ife2)
        self.assertTrue(ife1 <= ife2)
        self.assertTrue(ife2 > ife1)
        self.assertTrue(ife2 >= ife1)
        self.assertTrue(ife2 != ife1)


class IfeGroupTest(TestCase):
    def test_knows_if_structured_if_has_structured_chain(self):
        val = IfeGroup(IfeChain(chain='A', internal=5),
                       IfeChain(chain='B', internal=0))
        self.assertTrue(val.is_structured)

    def test_knows_if_structured_if_has_no_structured_chain(self):
        val = IfeGroup(IfeChain(chain='A', internal=4),
                       IfeChain(chain='B', internal=0))
        self.assertFalse(val.is_structured)

    def test_knows_how_many_chains_in_it(self):
        val = IfeGroup(IfeChain(chain='A', internal=4),
                       IfeChain(chain='B', internal=0))
        self.assertEquals(2, len(val))

    def test_has_an_id(self):
        val = IfeGroup(IfeChain(pdb='0111', chain='A', internal=5, model=10),
                       IfeChain(pdb='0111', chain='C', internal=0, model=10))
        self.assertEquals('0111|10|A', val.id)

    def test_uses_structured_only_in_id(self):
        val = IfeGroup(IfeChain(pdb='0111', chain='A', internal=5, model=10),
                       IfeChain(pdb='0111', chain='C', internal=0, model=0),
                       IfeChain(pdb='0111', chain='D', internal=6, model=10))
        self.assertEquals('0111|10|D+0111|10|A', val.id)

    def test_if_no_structured_in_id_uses_first(self):
        val = IfeGroup(IfeChain(pdb='0111', chain='A', internal=2, model=1),
                       IfeChain(pdb='0111', chain='C', internal=0, model=2),
                       IfeChain(pdb='0111', chain='D', internal=2, model=3))
        self.assertEquals('0111|1|A', val.id)

    def test_duplicate_additions_do_nothing(self):
        val = IfeGroup(IfeChain(pdb='0111', chain='A', internal=4, model=4),
                       IfeChain(pdb='0111', chain='C', internal=0, model=4))
        val.add(IfeChain(pdb='0111', chain='A', internal=4, model=4))
        self.assertEquals('0111|4|A', val.id)

    def test_dispatches_length_to_integral(self):
        val = IfeGroup(IfeChain(pdb='0111', chain='A', internal=4, length=5))
        self.assertEquals(5, val.length)

    def test_dispatches_full_length_to_integral(self):
        val = IfeGroup(IfeChain(pdb='0111', chain='A', internal=4,
                                full_length=5))
        self.assertEquals(5, val.full_length)

    def test_dispatches_internal_to_integral(self):
        val = IfeGroup(IfeChain(pdb='0111', chain='A', internal=4,
                                full_length=5))
        self.assertEquals(4, val.internal)

    def test_dispatches_pdb_to_integral(self):
        val = IfeGroup(IfeChain(pdb='0111', chain='A', internal=4,
                                full_length=5))
        self.assertEquals('0111', val.pdb)

    def test_dispatches_bps_to_integral(self):
        val = IfeGroup(IfeChain(pdb='0111', chain='A', internal=100,
                                full_length=5, bps=120),
                       IfeChain(pdb='0111', chain='B', internal=10,
                                full_length=5, bps=12))
        self.assertEquals(120, val.bps)

    def test_dispatches_model_to_integral(self):
        val = IfeGroup(IfeChain(pdb='0111', chain='A', internal=100,
                                full_length=5, model=2),
                       IfeChain(pdb='0111', chain='B', internal=10,
                                full_length=5, model=10))
        self.assertEquals(2, val.model)

    def test_getattr_complains_if_missing_key(self):
        val = IfeGroup(IfeChain(pdb='0111', chain='A', internal=4,
                                full_length=5))
        self.assertRaises(AttributeError, lambda: val.bob)

    def test_uses_most_bp_chain_as_integral(self):
        val = IfeGroup(IfeChain(pdb='0111', chain='A', internal=4, model=1),
                       IfeChain(pdb='0111', chain='C', internal=0, model=1))
        self.assertEquals('0111|1|A', val.integral.id)

    def test_uses_length_as_tiebreak(self):
        val = IfeGroup(IfeChain(pdb='0111', chain='A', length=10, internal=50,
                                model=10),
                       IfeChain(pdb='0111', chain='C', length=50, internal=50,
                                model=10))
        self.assertEquals('0111|10|C', val.integral.id)

    def test_uses_name_to_tiebreak_integral_chain(self):
        val = IfeGroup(IfeChain(pdb='0111', chain='C', length=10, internal=0,
                                model=3),
                       IfeChain(pdb='0111', chain='A', length=10, internal=0,
                                model=3))
        self.assertEquals('0111|3|A', val.integral.id)


class InfoLoadingTest(StageTest):
    loader_class = IfeLoader

    def test_loads_all_chains(self):
        val, _ = self.loader("4V4Q")
        self.assertEquals(6, len(val))

    def test_loads_interactions(self):
        _, val = self.loader("4V4Q")
        self.assertTrue(val)

    def test_loads_database_id(self):
        val = self.loader.load("4V4Q", "AA")
        assert isinstance(val.db_id, int)

    def test_loads_internal_cww(self):
        val = self.loader.load("4V4Q", "AA")
        self.assertEquals(472, val.internal)

    def test_loads_resolved_length(self):
        val = self.loader.load("4V4Q", "AA")
        self.assertEquals(1530, val.length)

    def test_loads_experimental_length(self):
        val = self.loader.load("4V4Q", "AA")
        self.assertEquals(1542, val.full_length)

    def test_has_an_id(self):
        val = self.loader.load("4V4Q", "AA", model=1)
        self.assertEquals("4V4Q|1|AA", val.id)

    def test_marks_yeast_5_8S_structured(self):
        val = self.loader.load("4V7R", "D1")
        self.assertTrue(val.is_structured)

    def test_can_load_correct_counts_with_sym_ops(self):
        val = self.loader.load('2QQP', 'R')
        self.assertEquals(val.length, 4)

    def test_can_load_bps_with_sym_ops(self):
        val = self.loader.load('1MDG', 'A')
        self.assertEquals(val.bps, 0)

    def test_can_load_correct_cww_with_sym_ops(self):
        val = self.loader.load('4PMI', 'A')
        self.assertEquals(val.length, 40)
        self.assertEquals(val.internal, 15)
        self.assertEquals(val.bps, 17)

    def test_can_load_chains_with_viral_sym_ops(self):
        ifes, inters = self.loader('1A34')
        names = [ife.chain for ife in ifes]
        bps = [ife.bps for ife in ifes]
        internal = [ife.internal for ife in ifes]
        self.assertEquals(['B', 'C'], names)
        self.assertEquals([0, 0], internal)
        self.assertEquals([0, 0], bps)
        ans = {'B': {'B': 0, 'C': 9}, 'C': {'B': 9, 'C': 0}}
        self.assertEquals(ans, inters)

    def test_loads_correct_length_4V9Q_BV(self):
        ife = self.loader.load('4V9Q', 'BV')
        self.assertEquals(77, ife.length)


class CrossChainInteractionsTest(StageTest):
    loader_class = IfeLoader

    def test_loads_interactions_symmetrically(self):
        ifes = [IfeChain(pdb='2MKN', chain='C'),
                IfeChain(pdb='2MKN', chain='B')]
        val = self.loader.cross_chain_interactions(ifes)
        ans = {'C': {'C': 0, 'B': 19}, 'B': {'B': 0, 'C': 19}}
        self.assertEquals(ans, val)

    def test_can_load_yeast_ribo_interactions(self):
        ifes = [IfeChain(pdb='4V7R', chain='D1'),
                IfeChain(pdb='4V7R', chain='D3')]
        val = self.loader.cross_chain_interactions(ifes)
        ans = {'D1': {'D1': 0, 'D3': 16}, 'D3': {'D1': 16, 'D3': 0}}
        self.assertEquals(ans, val)

    def test_can_load_correct_with_sym_ops(self):
        ifes = [IfeChain(pdb='1A34', chain='B'),
                IfeChain('1A34', chain='C')]
        val = self.loader.cross_chain_interactions(ifes, sym_op='P_1')
        ans = {'B': {'B': 0, 'C': 9}, 'C': {'B': 9, 'C': 0}}
        self.assertEquals(ans, val)


class SymOpTest(StageTest):
    loader_class = IfeLoader

    def test_if_one_picks_first(self):
        val = self.loader.sym_op("4V4Q")
        self.assertEquals(val, '1_555')

    def test_picks_a_sym_op(self):
        val = self.loader.sym_op('4PMI')
        self.assertEquals(val, '1_555')


class BestModelsTest(StageTest):
    loader_class = IfeLoader

    def test_with_one_model_picks_first(self):
        val = self.loader.best_model('4PMI', '1_555')
        self.assertEquals(val, 1)

    def test_it_picks_model_with_most_bps(self):
        val = self.loader.best_model('1E4P', '1_555')
        self.assertEquals(val, 8)
