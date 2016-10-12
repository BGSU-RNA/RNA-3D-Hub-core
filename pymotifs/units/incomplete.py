"""This is a loader to store data on which residues in a structure are
incomplete. Units can be incomplete in that they have atoms (or entire
residues) that have 0 occpancy. This is an issue that shows up in grouping
loops later.
"""

import itertools as it
import collections as coll

from pymotifs import core
from pymotifs import models as mod

from pymotifs.utils import row2dict

from pymotifs.units.info import Loader as InfoLoader

Entry = coll.namedtuple('Entry', ['model', 'chain', 'number', 'sequence',
                                  'alt_id', 'insertion'])


class Loader(core.SimpleLoader):
    allow_no_data = True
    dependencies = set([InfoLoader])

    def query(self, session, pdb):
        return session.query(mod.UnitIncomplete).\
            filter(mod.UnitIncomplete.pdb_id == pdb)

    def mapping(self, pdb):
        with self.session() as session:
            query = session.query(mod.UnitInfo.model,
                                  mod.UnitInfo.chain,
                                  mod.UnitInfo.number,
                                  mod.UnitInfo.unit.label('sequence'),
                                  mod.UnitInfo.alt_id,
                                  mod.UnitInfo.ins_code.label('insertion'),
                                  mod.UnitInfo.unit_id,
                                  ).\
                filter(mod.UnitInfo.pdb_id == pdb)
            mapping = coll.defaultdict(set)
            for result in query:
                r = row2dict(result)
                unit_id = r.pop('unit_id')
                key = Entry(**r)
                mapping[key].add(unit_id)
        return mapping

    def missing_keys(self, pdb):
        cif = self.cif(pdb)
        total = it.chain(getattr(cif, 'pdbx_unobs_or_zero_occ_residues', []),
                         getattr(cif, 'pdbx_unobs_or_zero_occ_atoms', []))
        missing = set()
        for row in total:
            model = int(row['PDB_model_num'])
            chain = row['auth_asym_id']
            num = int(row['auth_seq_id'])
            seq = row['auth_comp_id']
            ins = None
            if row['PDB_ins_code'] != '?':
                ins = row['PDB_ins_code']

            alt_id = None
            if 'label_alt_id' in row and row['label_alt_id'] != '?':
                alt_id = row['label_alt_id']

            key = Entry(model, chain, num, seq, alt_id, ins)
            missing.add(key)
        return missing

    def data(self, pdb, **kwargs):
        data = []
        mapping = self.mapping(pdb)
        for key in self.missing_keys(pdb):
            if key not in mapping:
                self.logger.info("Skipping unobserved %s", str(key))
            for uid in mapping[key]:
                data.append(mod.UnitIncomplete(unit_id=uid,
                                               pdb_id=pdb,
                                               is_incomplete=True))
        return data
