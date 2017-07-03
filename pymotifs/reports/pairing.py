"""
This is a module to produce a report about the
"""

from sqlalchemy.orm import aliased

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

    def positions(self, pdb, chain):
        with self.session() as session:
            esum = mod.ExpSeqUnitMapping
            esp = mod.ExpSeqPosition
            escm = mod.ExpSeqChainMapping
            query = session.query(
                esum.unit_id,
                esp.index,
                esp.unit,
            ).join(esp, esp.exp_seq_position_id == esum.exp_seq_position_id).\
                join(escm,
                     escm.exp_seq_chain_mapping_id == esum.exp_seq_chain_mapping_id).\
                join(mod.ChainInfo, mod.ChainInfo.chain_id == escm.chain_id).\
                filter(mod.ChainInfo.pdb_id == pdb).\
                filter(mod.ChainInfo.chain_name == chain)

            positions = []
            for result in query:
                entry = row2dict(result)
                entry['observed'] = int(result.unit_id is not None)
                positions.append(entry)
            return positions

    def interactions(self, pdb_id, chain, positions):
        mapping = {position['unit_id']: position for position in positions}
        with self.session() as session:
            uid1 = aliased(mod.UnitInfo)
            uid2 = aliased(mod.UnitInfo)
            query = session.query(mod.UnitPairsInteractions).\
                join(uid1,
                     uid1.unit_id == mod.UnitPairsInteractions.unit_id_1).\
                join(uid2,
                     uid2.unit_id == mod.UnitPairsInteractions.unit_id_2).\
                filter(mod.UnitPairsInteractions.f_lwbp == 'cWW').\
                filter(uid1.sym_op == uid2.sym_op)
            query = self.__limit_units__(query, uid1, pdb_id, chain)
            query = self.__limit_units__(query, uid2, pdb_id, chain)

            print(query)
            print(pdb_id)
            print(chain)
            print(query.count())
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
                filter(uid.sym_op.in_(['1_555', 'P_1']))
                # filter(uid.alt_id.in_([None, 'A'])).\

    def data(self, chain_specs, **kwargs):
        pdb, chain = chain_specs[0].split('.')
        positions = self.positions(pdb, chain)
        interactions = self.interactions(pdb, chain, positions)
        print(interactions)
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
