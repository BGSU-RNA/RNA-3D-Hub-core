"""
This is a module to produce a report about the
"""

from sqlachemly.orm import aliased

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import row2dict


class Reporter(core.Reporter):
    headers = [
        'unit_1',
        'index_1',
        'unit_id_1',
        'unit_2',
        'index_2',
        'unit_id_2',
        'observed',
    ]

    def exp_seq(self, pdb, chain):
        with self.session() as session:
            return session.query(mod.ExpSeqPdbMapping).\
                filter(pdb == pdb).\
                filter(chain == chain).\
                one().\
                exp_seq_id

    def load_positions(self, pdb, chain):
        exp_seq = self.exp_seq(pdb, chain)
        with self.session() as session:
            esum = mod.ExpSeqUnitMapping
            esp = mod.ExpSeqPosition
            query = session.query(
                esum.unit_id,
                esp.index,
                esp.unit,
            ).join(esp, esp.exp_seq_position_id == esum.exp_seq_position_id).\
                filter(esp.exp_seq_seq_id == exp_seq)

            positions = []
            for result in query:
                entry = row2dict(result)
                entry['observed'] = result['unit_id'] is not None
                positions.append(entry)
            return positions

    def load_interactions(self, pdb_id, chain, positions):
        mapping = {position['unit_id']: position for position in positions}
        with self.session() as session:
            uid1 = aliased(mod.UnitInfo)
            uid2 = aliased(mod.UnitInfo)
            query = session.query(mod.UnitPairsInteractions).\
                join(uid1,
                     uid1.unit_id == mod.UnitPairsInteractions.unit_id_1).\
                join(uid2,
                     uid2.unit_id == mod.UnitPairsInteractions.unit_id_2).\
                filter_by(f_lwbp='cWW').\
                filter(uid1.chain == uid2.chain).\
                filter(uid1.sym_op == uid2.sym_op).\
                filter(uid1.model == uid2.model)
            query = self.__limit_units__(query, uid1, pdb_id, chain)
            query = self.__limit_units__(query, uid2, pdb_id, chain)

            interactions = {}
            for result in query:
                unit = mapping[result.unit_id_1]['unit_id']
                interactions[unit] = mapping[result.unit_id_2]
            return interactions

    def __limit_units__(self, query, uid, pdb, chain):
        return query.\
                filter(uid.pdb_id == pdb).\
                filter(uid.chain == chain).\
                filter(uid.model == 1).\
                filter(uid.alt_id.in_([None, 'A'])).\
                filter(uid.sym_op.in_(['1_555', 'P_1']))

    def data(self, chain_spec):
        pdb, chain = chain_spec.split('.')
        positions = self.positions(pdb, chain)
        interactions = self.interactions(positions)
        for position in positions:
            base = {
                'unit_id_1': position['unit_id'],
                'unit_1': position['unit'],
                'index_1': position['index'],
                'observed': position['observed'],
            }
            other = interactions.get(position['unit_id'], {})
            second = {k + '_2': v for k, v in other.items() if k != 'observed'}
            base.update(second)
            yield base
