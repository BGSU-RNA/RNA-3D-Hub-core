import numpy as np

import pymotifs.core as core

from pymotifs.models import UnitInfo
from pymotifs.models import UnitPairsDistances
from pymotifs.units.info import Loader as InfoLoader


class Loader(core.SimpleLoader):
    max_distance = 10.0
    dependencies = set([InfoLoader])

    def query(self, session, pdb):
        return session.query(UnitPairsDistances).\
            join(UnitInfo, UnitInfo.id == UnitPairsDistances.unit1_id).\
            filter(UnitInfo.pdb_id == pdb)

    def center(self, residue):
        if 'base' in residue.centers:
            return residue.centers['base']
        if 'aa_backbone' in residue.centers:
            return residue.centers['backbone']
        if residue.sequence == 'HOH':
            return None
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
                if residue1 == residue2:
                    next

                distance = self.distance(residue1, residue2)
                if distance is not None and distance < self.max_distance:
                    yield UnitPairsDistances(unit1_id=residue1.unit_id(),
                                             unit2_id=residue2.unit_id(),
                                             distance=distance)
