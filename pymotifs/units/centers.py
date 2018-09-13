"""Load all unit centers into the database.
"""

import pymotifs.core as core

from pymotifs import models as mod

from pymotifs.units.info import Loader as InfoLoader


class Loader(core.SimpleLoader):
    """A class to load all base centers into the database.
    """

    dependencies = set([InfoLoader])
    allow_no_data = True

    def query(self, session, pdb):
        return session.query(mod.UnitCenters).\
            filter(mod.UnitCenters.pdb_id == pdb)

    def data(self, pdb, **kwargs):
        structure = self.structure(pdb)
        for residue in structure.residues():
            for name in residue.centers.definitions():
                center = residue.centers[name]
                if len(center) == 3:
                    yield mod.UnitCenters(unit_id=residue.unit_id(),
                                          name=name,
                                          pdb_id=pdb,
                                          x=float(center[0]),
                                          y=float(center[1]),
                                          z=float(center[2]))
