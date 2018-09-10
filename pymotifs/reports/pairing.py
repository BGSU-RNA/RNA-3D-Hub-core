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
        'unit_2',
        'index_2',
        'observed',
    ]

    def exp_seq(self, pdb_id, chain):
        with self.session() as session:
            return session.query(mod.ExpSeqPdb).\
                filter_by(pdb_id=pdb_id, chain_name=chain).\
                one().\
                exp_seq_id

    def positions(self, pdb, chain):
        exp_seq = self.exp_seq(pdb, chain)
        with self.session() as session:
            esum = mod.ExpSeqUnitMapping
            esp = mod.ExpSeqPosition
            escm = mod.ExpSeqChainMapping
            ci = mod.ChainInfo
            query = session.query(
                esum.unit_id,
                esp.index,
                esp.unit,
            ).join(esp, esp.exp_seq_position_id == esum.exp_seq_position_id).\
                join(escm, escm.exp_seq_chain_mapping_id == esum.exp_seq_chain_mapping_id).\
                join(ci, ci.chain_id == escm.chain_id).\
                filter(ci.pdb_id == pdb).\
                filter(ci.chain_name == chain)

            if not query.count():
                raise core.InvalidState("Could not load positions for %s|1|%s" % (pdb_id, chain))

            positions = []
            for result in query:
                entry = row2dict(result)
                entry['observed'] = int(result.unit_id is not None)
                entry['index'] = entry['index'] + 1
                positions.append(entry)
            return positions

    def interactions(self, pdb_id, chain, positions, remove_pseudoknots=False):
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

            if remove_pseudoknots:
                query = query.filter(mod.UnitPairsInteractions.f_crossing < 4)

            if not query.count():
                raise core.InvalidState("Could not load interactions for %s|1|%s" % (pdb_id, chain))

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

    def data(self, chain_spec, remove_pseudoknots, **kwargs):
        pdb, chain = chain_spec
        positions = self.positions(pdb, chain)
        interactions = self.interactions(pdb, chain, positions,
                                         remove_pseudoknots=remove_pseudoknots)
        first_only = set(['observed', 'unit_id'])
        for position in positions:
            base = {
                'unit_1': position['unit'],
                'index_1': position['index'],
                'observed': position['observed'],
            }
            other = interactions.get(position['unit_id'], {})
            second = {k + '_2': v for k, v in other.items() if k not in first_only}
            base.update(second)
            yield base
