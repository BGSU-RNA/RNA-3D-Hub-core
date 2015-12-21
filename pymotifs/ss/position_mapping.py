from pymotifs import core
from pymotifs import utils as ut
from pymotifs import models as mod

from pymotifs.utils.alignment import align

from pymotifs.ss.exp_seq_mapping import Loader as SsMappingLoader
from pymotifs.ss.positions import Loader as SsPositionLoader
from pymotifs.exp_seq.positions import Loader as ExpPositionLoader
from pymotifs.exp_seq.mapping import Loader as ExpSeqMappingLoader
from pymotifs.correspondence.positions import Loader as CorrPositionLoader


class Loader(core.SimpleLoader):
    dependencies = set([SsMappingLoader, SsPositionLoader, ExpPositionLoader,
                        ExpSeqMappingLoader])
    table = mod.SsExpSeqPositionMapping

    def to_process(self, pdbs, **kwargs):
        ss_map = mod.SsExpSeqMapping
        exp_pdbs = mod.ExpSeqPdb
        seqs = set()
        for chunk in ut.grouper(1000, pdbs):
            with self.session() as session:
                query = session.query(ss_map.ss_exp_seq_mapping_id).\
                    join(exp_pdbs, ss_map.exp_seq_id == exp_pdbs.exp_seq_id).\
                    filter(exp_pdbs.pdb_id.in_(chunk)).\
                    distinct()
                seqs.update(r.ss_exp_seq_mapping_id for r in query)
        return sorted(seqs)

    def query(self, session, map_id):
        return session.query(self.table).\
            filter_by(ss_exp_seq_mapping_id=map_id)

    def pair(self, map_id):
        with self.session() as session:
            result = session.query(mod.SsExpSeqMapping).\
                filter_by(ss_exp_seq_mapping_id=map_id).\
                one()
            return (result.exp_seq_id, result.ss_id)

    def exp_info(self, exp_id):
        loader = CorrPositionLoader(self.config, self.session)
        return loader.sequence(exp_id)

    def ss_info(self, ss_id):
        with self.session() as session:
            query = session.query(mod.SsPositions).\
                filter_by(ss_id=ss_id).\
                order_by(mod.SsPositions.index)

            ids = []
            sequence = []
            for result in query:
                ids.append(result.ss_position_id)
                sequence.append(result.unit)
        return {'ids': ids, 'sequence': ''.join(sequence)}

    def align(self, exp_info, ss_info):
        alignment = align([exp_info, ss_info])
        return [{'ss': ss_id, 'exp': exp_id} for exp_id, ss_id in alignment]

    def data(self, map_id, **kwargs):
        exp_id, ss_id = self.pair(map_id)
        exp = self.exp_info(exp_id)
        ss = self.ss_info(ss_id)
        alignment = self.align(exp, ss)
        data = []
        for position in alignment:
            data.append({
                'ss_position_id': position['ss'],
                'exp_seq_position_id': position['exp'],
                'ss_exp_seq_mapping_id': map_id
            })
        return data
