import pymotifs.core as core

from pymotifs.models import UnitInfo
from pymotifs.models import UnitPairsDistances


class Loader(core.Loader):
    max_distance = 10.0

    def transform(self, pdb):
        return [self.cif(pdb)]

    def has_data(self, structure):
        with self.session() as session:
            query = session.query(UnitPairsDistances).\
                join(UnitInfo, UnitInfo.id == UnitPairsDistances.unit1_id).\
                filter(UnitInfo.pdb == structure.pdb)
            return bool(query.count())

    def remove(self, structure):
        with self.session() as session:
            session.query(UnitPairsDistances).\
                join(UnitInfo, UnitInfo.id == UnitPairsDistances.unit1_id).\
                filter(UnitInfo.pdb == structure.pdb).\
                delete()

    def data(self, structure):
        for residue1 in structure.residues():
            for residue2 in structure.residues():
                if residue1 == residue2 or residue1 < residue2:
                    next
                distance = residue1.distance(residue2)
                if distance < self.max_distance:
                    yield UnitPairsDistances(unit1_id=residue1.unit_id(),
                                             unit2_id=residue2.unit_id(),
                                             distance=distance)
