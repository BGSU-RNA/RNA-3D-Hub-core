import pymotifs.core as core

from pymotifs.models import UnitInfo
from pymotifs.models import UnitCoordinates
from pymotifs.units.info import Loader as InfoLoader


class Loader(core.SimpleLoader):
    dependencies = set([InfoLoader])

    def query(self, session, pdb):
        return session.query(UnitCoordinates).\
            join(UnitInfo, UnitInfo.id == UnitCoordinates.id).\
            filter(UnitInfo.pdb_id == pdb)

    def data(self, pdb):
        pass
