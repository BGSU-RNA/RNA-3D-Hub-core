"""Load all unit centers into the database.
"""

import pymotifs.core as core

from pymotifs.models import UnitInfo
from pymotifs.models import UnitCenters

from pymotifs.units.info import Loader as InfoLoader


class Loader(core.SimpleLoader):
    """A class to load all base centers into the database.
    """

    dependencies = set([InfoLoader])
    allow_no_data = False

    def query(self, session, pdb):
        return session.query(UnitCenters).\
            join(UnitInfo, UnitInfo.unit_id == UnitCenters.unit_id).\
            filter(UnitInfo.pdb_id == pdb)

    def data(self, pdb, **kwargs):
        structure = self.structure(pdb)
        for residue in structure.residues():
            for name in residue.centers.definitions():
                center = residue.centers[name]
                if len(center) == 3:
                    yield UnitCenters(unit_id=residue.unit_id(),
                                      name=name,
                                      pdb_id=pdb,
                                      x=center[0],
                                      y=center[1],
                                      z=center[2])
