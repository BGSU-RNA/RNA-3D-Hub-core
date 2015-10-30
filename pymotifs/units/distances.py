import numpy as np

import pymotifs.core as core

from pymotifs.models import UnitInfo
from pymotifs.models import UnitPairsDistances
from pymotifs.units.info import Loader as InfoLoader


class Loader(core.Loader):
    max_distance = 10.0
    dependencies = set([InfoLoader])

    def remove(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(UnitPairsDistances.unit_id_1.label('unit_id')).\
                join(UnitInfo, UnitInfo.unit_id == UnitPairsDistances.unit_id_1).\
                filter(UnitInfo.pdb_id == pdb)
            ids = [result.unit_id for result in query]

        with self.session() as session:
            query = session.query(UnitPairsDistances).\
                filter(UnitPairsDistances.unit_id_1.in_(ids)).\
                delete(synchronize_session=True)

    def has_data(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(UnitPairsDistances).\
                join(UnitInfo, UnitInfo.unit_id == UnitPairsDistances.unit_id_1).\
                filter(UnitInfo.pdb_id == pdb)

            return bool(query.count())

    def center(self, residue):
        if residue.sequence == 'HOH':
            return None

        names = ['base', 'aa_backbone']
        for name in names:
            if name in residue.centers:
                return residue.centers[name]

        return np.mean(residue.coordinates(), axis=0)

    def distance(self, residue1, residue2):
        center1 = self.center(residue1)
        center2 = self.center(residue2)

        if center1 is not None and center2 is not None:
            return np.linalg.norm(center1 - center2)
        return None

    def data(self, pdb, **kwargs):
        structure = self.structure(pdb)
        for residue1 in structure.residues():
            for residue2 in structure.residues():
                if residue1.unit_id() == residue2.unit_id():
                    next

                distance = self.distance(residue1, residue2)
                if distance is not None and distance < self.max_distance:
                    yield UnitPairsDistances(unit_id_1=residue1.unit_id(),
                                             unit_id_2=residue2.unit_id(),
                                             distance=distance)
