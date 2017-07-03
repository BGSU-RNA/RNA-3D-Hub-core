"""
This is a module to produce a report about the
"""

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import row2dict


class Reporter(core.Reporter):
    headers = [
        'unit',
        'interacting_index',
        'observed',
    ]

# select
# 	UI1.unit_id,
# 	ESP1.index,
# 	ESP1.unit,
# 	UI2.unit_id,
# 	ESP2.index,
# 	ESP2.unit,
# 	UPI.f_lwbp
# from exp_seq_unit_mapping as ESUM1
# join exp_seq_position as ESP1
# on
# 	ESP1.exp_seq_position_id = ESUM1.exp_seq_position_id
# left join unit_info as UI1
# on
# 	UI1.unit_id = ESUM1.unit_id
# left join unit_pairs_interactions as UPI
# on
# 	UPI.unit_id_1 = UI1.unit_id
# left join unit_info as UI2
# on
# 	UI2.unit_id = UPI.unit_id_2
# left join exp_seq_unit_mapping as ESUM2
# on
# 	ESUM2.unit_id = UI2.unit_id
# left join exp_seq_position as ESP2
# on
# 	ESP2.exp_seq_position_id = ESUM2.exp_seq_position_id
# where
# 	UI1.pdb_id = '5AOX'
# 	and UI1.chain = 'C'
# 	and UPI.f_lwbp = 'cWW'
# order by ESP1.index
# ;

    def load_interactions(self, pdb, chain, model=1):
        with self.session() as session:
            ui1 = aliased(mod.UnitInfo)
            ui2 = aliased(mod.UnitInfo)
            esp1 = aliased(mod.ExpSeqPosition)
            esp2 = aliased(mod.ExpSeqPosition)
            esum1 = aliased(mod.ExpSeqUnitMapping)
            esum2 = aliased(mod.ExpSeqUnitMapping)
            upi = mod.UnitPairsInteractions

            query = session.query(
                esum1.unit_id,
                esp1.index,
                esp1.unit,
                ui2.unit_id,
                esp2.index,
                esp2.unit,
            ).\
                outerjoin(esp1,
                          esp.exp_seq_position_id == esum1.exp_seq_position_id).\
                outerjoin(ui1, ui1.unit_id == esum1.unit_id).\
                outerjoin(upi, upi.unit_id_1 == ui1.unit_id).\
                outerjoin(ui2, ui2.unit_id == upi.unit_id_2).\
                outerjoin(esum2,
                          esum2.unit_id == ui2.unit_id).\
                outerjoin(esp2,
                          esp2.exp_seq_position_id == esum2.exp_seq_position_id).\
                filter(upi.f_lwbp == 'cWW').\
                filter(ui1.chain == ui2.chain).\
                filter(ui1.sym_op == ui2.sym_op).\
                filter(ui1.pdb == pdb).\
                filter(ui1.chain == chain).\
                filter(ui1.sym_op.in_(('1_555', 'P_1'))).\
                order_by(esp1.index)

            return [row2dict(r) for r in query]

    def data(self, chain_spec):
        pdb, chain = chain_spec.split('.')
        return self.load_interactions(pdb, chain, model=1)
