import numpy as np

import pymotifs.core as core

from pymotifs.models import UnitInfo
from pymotifs.models import UnitPairsDistances


class Loader(core.SimpleLoader):
    max_distance = 10.0

    def transform(self, pdb):
        return [self.cif(pdb)]

    def query(self, session, structure):
        return session.query(UnitPairsDistances).\
            join(UnitInfo, UnitInfo.id == UnitPairsDistances.unit1_id).\
            filter(UnitInfo.pdb == structure.pdb)

    def distance(self, residue1, residue2):
        return np.linalg.norm(residue1.centers['center'],
                              residue2.centers['center'])

    def data(self, structure):
        for residue1 in structure.residues():
            for residue2 in structure.residues():
                if residue1 == residue2 or residue1 > residue2:
                    next

                distance = self.distance(residue1, residue2)
                if distance < self.max_distance:
                    yield UnitPairsDistances(unit1_id=residue1.unit_id(),
                                             unit2_id=residue2.unit_id(),
                                             distance=distance)
