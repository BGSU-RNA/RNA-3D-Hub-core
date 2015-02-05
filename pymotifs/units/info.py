import pymotifs.core as core

from Bio.Alphabet import ThreeLetterProtein

from pymotifs.models import UnitInfo

AA = [seq.upper() for seq in ThreeLetterProtein().letters]


class Loader(core.SimpleLoader):
    update_gap = False

    def query(self, session, pdb):
        return session(UnitInfo).filter_by(pdb_id=pdb)

    def transform(self, pdb):
        return self.cif(pdb).structure()

    def type(self, unit):
        seq = unit.sequence.upper()
        if seq in ['A', 'C', 'G', 'U']:
            return 'rna'
        if seq == 'HOH':
            return 'water'
        if seq in AA:
            return 'aa'
        if seq in ['DA', 'DC', 'DG', 'DT']:
            return 'dna'
        return None

    def as_unit(self, nt):
        return UnitInfo(id=nt.unit_id(),
                        pdb_id=nt.pdb,
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
