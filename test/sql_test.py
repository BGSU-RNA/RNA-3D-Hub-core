from unittest import TestCase

from sqlalchemy import Column, Integer, String, ForeignKey
from sqlalchemy.orm import aliased

from pymotifs import models as mod
from pymotifs.core import Session as Wrapper

from test import Session


class SimpleTypes(mod.Base):
    __tablename__ = '___testing_type'
    __table_args__ = {'useexisting': True}

    id = Column(String(10), primary_key=True)
    type_name = Column(String(50))


class Simple(mod.Base):
    __tablename__ = '___testing_simple'
    __table_args__ = {'useexisting': True}

    id = Column(Integer, primary_key=True)
    name = Column(String(50))
    pdb_id = Column(String(4))
    type_id = Column(String(10), ForeignKey('___testing_type.id'))


class DeleteWithJoinsTest(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tables = [SimpleTypes.__table__, Simple.__table__]
        mod.Base.metadata.create_all(tables=cls.tables)

    @classmethod
    def tearDownClass(cls):
        mod.Base.metadata.drop_all(tables=cls.tables)

    def setUp(self):
        self.session = Wrapper(Session)
        self.fixtures()

    def tearDown(self):
        with self.session() as session:
            session.execute('delete from ___testing_simple')
            session.execute('delete from ___testing_type')

    def fixtures(self):
        with self.session() as session:
            session.add_all([
                SimpleTypes(id='Person', type_name='For people'),
                SimpleTypes(id='PDB', type_name='ID of a pdb'),
            ])

        with self.session() as session:
            session.add_all([
                Simple(id=3, name='BOB', pdb_id=None, type_id='Person'),
                Simple(id=4, name='Euk SSU', pdb_id='4V88', type_id='PDB'),
                Simple(id=5, name='Tt SSU', pdb_id='4V8G', type_id='PDB'),
                Simple(id=6, name='Tt SSU', pdb_id='4V9Q', type_id='PDB'),
            ])

    def assertHasIds(self, table, ids):
        with self.session() as session:
            query = session.query(table.id)
            val = [result.id for result in query]
        self.assertEquals(sorted(val), sorted(ids))

    def test_can_use_simple_deletes(self):
        with self.session() as session:
            session.query(Simple).delete(synchronize_session=False)
        self.assertHasIds(Simple, [])
        self.assertHasIds(SimpleTypes, ['Person', 'PDB'])

    def test_can_use_delete_with_where(self):
        with self.session() as session:
            session.query(Simple).\
                filter(Simple.name == 'Tt SSU').\
                delete(synchronize_session=False)
        self.assertHasIds(Simple, [3, 4])
        self.assertHasIds(SimpleTypes, ['Person', 'PDB'])

    def test_it_can_do_nothing(self):
        with self.session() as session:
            session.query(Simple).\
                filter(Simple.name == 'NOTHING').\
                delete(synchronize_session=False)
        self.assertHasIds(Simple, [3, 4, 5, 6])
        self.assertHasIds(SimpleTypes, ['Person', 'PDB'])

    def test_can_use_delete_with_aliases_and_where(self):
        with self.session() as session:
            simple = aliased(Simple)
            session.query(simple).\
                filter(simple.name == 'Tt SSU').\
                delete(synchronize_session=False)
        self.assertHasIds(Simple, [3, 4])
        self.assertHasIds(SimpleTypes, ['Person', 'PDB'])

    # def test_can_use_delete_with_join(self):
    #     with self.session() as session:
    #         simple = aliased(Simple)
    #         session.query(simple).\
    #             join(SimpleTypes, SimpleTypes.id == simple.type_id).\
    #             filter(SimpleTypes.type_name == 'For people').\
    #             delete(synchronize_session=False)
    #     self.assertCount(Simple, [4, 5, 6])
    #     self.assertHasIds(SimpleTypes, ['Person', 'PDB'])

    # def test_can_use_complex_join(self):
    #     with self.session() as session:
    #         session.query(Simple).\
    #             join(SimpleTypes, SimpleTypes.type_name == 'ID of a pdb').\
    #             join(mod.PdbInfo, mod.PdbInfo.pdb_id == Simple.pdb_id).\
    #             join(mod.ChainInfo, mod.ChainInfo.pdb_id == mod.PdbInfo.pdb_id).\
    #             join(mod.UnitInfo,
    #                  (mod.UnitInfo.chain == mod.ChainInfo.chain_name) &
    #                  (mod.UnitInfo.pdb_id == mod.PdbInfo.pdb_id)).\
    #             filter(mod.ChainInfo.chain_name.in_(['AA', '-G'])).\
    #             filter(mod.UnitInfo.unit_type_id != 'rna').\
    #             delete(synchronize_session=False)
    #     self.assertCount(Simple, [3, 5, 6])
    #     self.assertHasIds(SimpleTypes, ['Person', 'PDB'])
