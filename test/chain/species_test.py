import pytest

from test import StageTest

from pymotifs.chains.species import Loader
from pymotifs.models import ChainInfo


class QueryingTest(StageTest):
    loader_class = Loader

    def test_knows_has_no_data(self):
        self.assertFalse(self.loader.has_data('0GID'))

    def test_knows_has_data(self):
        self.assertTrue(self.loader.has_data('124D'))

    @pytest.mark.skip()
    def test_can_remove_data(self):
        pass


class GettingDataTest(StageTest):
    loader_class = Loader

    def test_can_get_data_for_pdb(self):
        data = self.loader.data('1GID')
        val = [d['species_id'] for d in data]
        ans = [32630, 32630]
        self.assertEquals(ans, val)

    def test_assigns_none_when_unknown(self):
        data = self.loader.data('1ET4')
        val = [d['species_id'] for d in data]
        ans = [None] * 5
        self.assertEquals(ans, val)

    def test_assigns_to_species(self):
        data = self.loader.data('4V6M')
        val = [e['species_id'] for e in data]
        ans = [562] * 5
        self.assertEquals(ans, val)

    def test_can_find_species_when_name_and_id_differ(self):
        data = self.loader.data('4V9Q')
        data.sort(key=lambda c: c['chain_id'])
        val = [entry['species_id'] for entry in data]
        #      AA,  AB,  BA,  BV,  BW,  BX,    CA,  CB,  DA,  DV,  DW,  DX
        ans = [274, 274, 274, 512, 512, 32630, 274, 274, 274, 512, 512, 32630]
        self.assertEquals(ans, val)


class ValidatingTest(StageTest):
    loader_class = Loader

    def id_of(self, pdb, name):
        with self.loader.session() as session:
            return session.query(ChainInfo).\
                filter_by(pdb_id=pdb, chain_name=name).\
                one().\
                chain_id

    def test_will_try_to_correct_bad_ecoli(self):
        chain_id = self.id_of('4V9Q', 'BV')
        val = self.loader.validated(chain_id, None)
        self.assertEquals(562, val)

    def test_will_try_to_correct_bad_512(self):
        chain_id = self.id_of('4V9Q', 'BV')
        val = self.loader.validated(chain_id, 512)
        self.assertEquals(562, val)

    def test_will_return_good_name(self):
        chain_id = self.id_of('1FJG', 'A')
        val = self.loader.validated(chain_id, 274)
        self.assertEquals(274, val)
