import pytest

from test import StageTest

from pymotifs.units.info import Loader

from fr3d.cif.reader import Cif
from fr3d.data import Component


class DetectingComponentTypeTest(StageTest):
    loader_class = Loader

    def test_knows_if_unit_is_rna(self):
        val = self.loader.type(Component([], sequence='A'))
        ans = 'rna'
        self.assertEqual(ans, val)

    def test_knows_if_unit_is_dna(self):
        val = self.loader.type(Component([], sequence='DT'))
        ans = 'dna'
        self.assertEqual(ans, val)

    def test_knows_if_unit_is_aa(self):
        val = self.loader.type(Component([], sequence='GLU'))
        ans = 'aa'
        self.assertEqual(ans, val)

    def test_knows_if_unit_is_water(self):
        val = self.loader.type(Component([], sequence='HOH'))
        ans = 'water'
        self.assertEqual(ans, val)

    def test_gives_none_for_unknown_unit(self):
        val = self.loader.type(Component([], sequence='GTP'))
        ans = None
        self.assertEqual(ans, val)


class CreatingUnitsTest(StageTest):
    loader_class = Loader

    @classmethod
    def setUpClass(cls):
        with open('test/files/cif/1GID.cif', 'rb') as raw:
            cls.structure = Cif(raw).structure()

    def setUp(self):
        super(CreatingUnitsTest, self).setUp()
        residues = self.__class__.structure.residues(number=103, chain='A')
        self.residues = list(residues)
        self.data = self.loader.as_unit(self.residues[0])

    def test_sets_the_id(self):
        val = self.data.unit_id
        ans = '1GID|1|A|G|103'
        self.assertEqual(ans, val)

    def test_sets_the_pdb_id(self):
        val = self.data.pdb_id
        ans = '1GID'
        self.assertEqual(ans, val)

    def test_sets_the_model(self):
        self.assertEquals(1, self.data.model)

    def test_sets_the_chain(self):
        val = self.data.chain
        ans = 'A'
        self.assertEqual(ans, val)

    def test_sets_the_unit(self):
        val = self.data.unit
        ans = 'G'
        self.assertEqual(ans, val)

    def test_sets_the_number(self):
        val = self.data.number
        ans = 103
        self.assertEqual(ans, val)

    def test_sets_the_alt_id(self):
        val = self.data.alt_id
        ans = None
        self.assertEqual(ans, val)

    def test_sets_the_ins_code(self):
        val = self.data.ins_code
        ans = None
        self.assertEqual(ans, val)

    def test_sets_the_sym_op(self):
        val = self.data.sym_op
        ans = '1_555'
        self.assertEqual(ans, val)

    @pytest.mark.xfail(reason='Not computing index yet')
    def test_sets_the_chain_index(self):
        assert self.data.chain_index == 0

    @pytest.mark.xfail(reason='Not computing index yet')
    def test_sets_the_global_index(self):
        assert self.data.global_index == 0

    def test_sets_the_type(self):
        val = self.data.unit_type_id
        ans = 'rna'
        self.assertEqual(ans, val)


class BuildingAllUnitsTest(StageTest):
    loader_class = Loader

    @classmethod
    def setUpClass(cls):
        with open('test/files/cif/1GID.cif', 'rb') as raw:
            cls.structure = Cif(raw).structure()

    def setUp(self):
        super(BuildingAllUnitsTest, self).setUp()
        self.structure = self.__class__.structure
        self.data = list(self.loader.data(self.structure))

    def test_loads_all_units(self):
        val = len(self.data)
        ans = 350
        self.assertEqual(ans, val)


class QueryingTest(StageTest):
    loader_class = Loader

    def test_creates_query_that_counts(self):
        with self.loader.session() as session:
            query = self.loader.query(session, '1GID')
            self.assertEquals(350, query.count())

    def test_knows_data_exists(self):
        self.assertTrue(self.loader.has_data('1GID'))

    def test_knows_if_no_data_exists(self):
        self.assertFalse(self.loader.has_data('0GID'))
