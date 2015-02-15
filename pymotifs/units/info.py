import pymotifs.core as core

from Bio.Alphabet import ThreeLetterProtein

from pymotifs.models import UnitInfo
from pymotifs.utils import units

AA = [seq.upper() for seq in ThreeLetterProtein().letters]


class Loader(core.SimpleLoader):
    update_gap = False

    def query(self, session, structure):
        return session.query(UnitInfo).filter_by(pdb_id=structure.pdb)

    def transform(self, pdb, **kwargs):
        return [self.cif(pdb).structure()]

    def type(self, unit):
        return units.component_type(unit)

    def as_unit(self, nt):
        return UnitInfo(id=nt.unit_id(),
                        pdb_id=nt.pdb,
                        model=nt.model,
                        chain=nt.chain,
                        unit=nt.sequence,
                        number=nt.number,
                        alt_id=getattr(nt, 'alt_id', None),
                        ins_code=nt.insertion_code,
                        sym_op=nt.symmetry,
                        chain_index=getattr(nt, 'chain_index', None),
                        unit_type_id=self.type(nt))

    def data(self, structure):
        return [self.as_unit(nt) for nt in structure.residues(polymeric=None)]
