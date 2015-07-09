import pymotifs.core as core

from pymotifs.models import UnitInfo
from pymotifs.models import UnitCenters


class Loader(core.SimpleLoader):

    def query(self, session, pdb):
        return session.query(UnitCenters).\
            join(UnitInfo, UnitInfo.id == UnitCenters.id).\
            filter(UnitInfo.pdb == pdb)

    def data(self, pdb, **kwargs):
        structure = self.structure(pdb)
        for residue in structure.residues():
            if residue.type == 'rna':
                center = residue.centers['base']
                yield UnitCenters(id=residue.unit_id(),
                                  name='base',
                                  x=center[0],
                                  y=center[1],
                                  z=center[2])
