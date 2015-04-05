import numpy as np

import pymotifs.core as core
from pymotifs.utils import units

from pymotifs.models import UnitInfo
from pymotifs.models import UnitPairsDistances
from pymotifs.units.info import Loader as InfoLoader


class Loader(core.SimpleLoader):
    max_distance = 10.0
    dependencies = set(InfoLoader)

    def query(self, session, pdb):
        return session.query(UnitPairsDistances).\
            join(UnitInfo, UnitInfo.id == UnitPairsDistances.unit1_id).\
            filter(UnitInfo.pdb_id == pdb)

    def center(self, residue):
        type = units.component_type(residue)
        if type == 'rna' or type == 'dna':
            return residue.centers['base']
        if type == 'aa':
            return residue.centers['backbone']
        if type == 'water' or type is None:
            return None
        return residue.center

    def distance(self, residue1, residue2):
        center1 = self.center(residue1)
        center2 = self.center(residue2)

        if center1 and center2:
            return np.linalg.norm(center1, center2)
        return None

    def data(self, pdb):
        structure = self.structure(pdb)
        for residue1 in structure.residues():
            for residue2 in structure.residues():
                if residue1 == residue2 or residue1 > residue2:
                    next

                distance = self.distance(residue1, residue2)
                if distance is not None and distance < self.max_distance:
                    yield UnitPairsDistances(unit1_id=residue1.unit_id(),
                                             unit2_id=residue2.unit_id(),
                                             distance=distance)
