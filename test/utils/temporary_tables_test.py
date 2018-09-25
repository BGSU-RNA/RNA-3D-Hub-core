import functools as ft

from test import StageTest

from pymotifs.models import TempPdbs
from pymotifs.utils.temporary_tables import Temporary


class TempTableTest(StageTest):
    loader_class = ft.partial(Temporary, TempPdbs)

    def test_can_build_a_temp_table(self):
        with self.loader() as session:
            query = session.query(TempPdbs.pdb_id).distinct()
            self.assertEquals([], [result.pdb_id for result in query])

    def test_removes_temp_table_afterwards(self):
        pdbs = ['4V4Q', '1FJG']
        with self.loader({'pdb_id': pdb} for pdb in pdbs) as session:
            pass

        try:
            val = []
            with self.loader.session() as session:
                query = session.query(TempPdbs.pdb_id).distinct()
                val = [result.pdb_id for result in query]
                self.fail("Query should fail from missing table")
        except Exception:
            pass
        finally:
            self.assertEquals([], val)

    def test_can_build_table_with_data(self):
        pdbs = ['1FJG', '4V4Q']
        with self.loader({'pdb_id': pdb} for pdb in pdbs) as session:
            query = session.query(TempPdbs.pdb_id).\
                distinct().\
                order_by(TempPdbs.pdb_id)
            val = [result.pdb_id for result in query]
        self.assertEquals(pdbs, val)
