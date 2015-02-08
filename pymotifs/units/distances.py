import numpy as np

import pymotifs.core as core

from pymotifs.models import UnitInfo
from pymotifs.models import UnitPairsDistances


class Loader(core.SimpleLoader):
    max_distance = 10.0

    def transform(self, pdb, **kwargs):
        return [self.structure(pdb)]

    def query(self, session, structure):
        return session.query(UnitPairsDistances).\
            join(UnitInfo, UnitInfo.id == UnitPairsDistances.unit1_id).\
            filter(UnitInfo.pdb == structure.pdb)

    def center(self, residue):
        if residue.type == 'rna':
            return residue.centers['base']
        if residue.type == 'aa':
            return residue.centers['backbone']
        if residue.type == 'dna':
            return residue.centers['base']
        return None

    def distance(self, residue1, residue2):
        center1 = self.center(residue1)
        center2 = self.center(residue2)

        if center1 and center2:
            return np.linalg.norm(center1, center2)
        return None

    def data(self, structure):
        for residue1 in structure.residues():
            for residue2 in structure.residues():
                if residue1 == residue2 or residue1 > residue2:
                    next

                distance = self.distance(residue1, residue2)
                if distance is not None and distance < self.max_distance:
                    yield UnitPairsDistances(unit1_id=residue1.unit_id(),
                                             unit2_id=residue2.unit_id(),
                                             distance=distance)
