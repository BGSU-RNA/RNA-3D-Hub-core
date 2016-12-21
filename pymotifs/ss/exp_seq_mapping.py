import os

from sqlalchemy.orm import aliased

from pymotifs import core
from pymotifs import utils as ut
from pymotifs import models as mod

from pymotifs.ss.info import Loader as SsInfoLoader
from pymotifs.exp_seq.chain_mapping import Loader as ExpSeqChainLoader
from pymotifs.ife.info import Loader as IfeLoader

from pymotifs.ss.helpers import md5


class Loader(core.SimpleLoader):
    dependencies = set([SsInfoLoader, IfeLoader, ExpSeqChainLoader])
    @property
    def table(self):
        return mod.SsExpSeqMapping
    mark = False

    def to_process(self, pdbs, **kwargs):
        seqs = set()
        for chunk in ut.grouper(100, pdbs):
            with self.session() as session:
                query = session.query(mod.ExpSeqPdb).\
                    filter(mod.ExpSeqPdb.pdb_id.in_(chunk)).\
                    distinct()

                if kwargs.get('chain'):
                    query = query.filter_by(chain_name=kwargs.get('chain'))

                seqs.update(result.exp_seq_id for result in query)
        return sorted(seqs)

    def query(self, session, exp_seq):
        return session.query(self.table).\
            filter_by(exp_seq_id=exp_seq)

    def ss_id(self, filename):
        os.chdir(self.config['locations']['base'])
        with self.session() as session:
            return session.query(mod.SsInfo).\
                filter_by(md5=md5(filename)).\
                one().\
                ss_id

    def mapped_id(self, exp_seq):
        if not self.has_data(exp_seq):
            return None

        with self.session() as session:
            return self.query(session, exp_seq).\
                one().\
                ss_id

    def nr_seqs(self, exp_seq):
        exp = aliased(mod.ExpSeqInfo)
        exp_map = aliased(mod.ExpSeqChainMapping)
        ife = aliased(mod.IfeChains)
        nr = aliased(mod.NrChains)
        mem = aliased(mod.NrChains)
        mem_ife = aliased(mod.IfeChains)
        mem_map = aliased(mod.ExpSeqChainMapping)
        with self.session() as session:
            query = session.query(mem_map.exp_seq_id).\
                join(mem_ife, mem_ife.chain_id == mem_map.chain_id).\
                join(mem, mem.ife_id == mem_ife.ife_id).\
                join(nr, nr.nr_class_id == mem.nr_class_id).\
                join(ife, ife.ife_id == nr.ife_id).\
                join(exp_map, exp_map.chain_id == ife.chain_id).\
                join(exp, exp.exp_seq_id == exp_map.exp_seq_id).\
                filter(exp.exp_seq_id == exp_seq).\
                distinct()
            return set(result.exp_seq_id for result in query)

    def matching_sequence(self, exp_seq):
        aligned = self.nr_seqs(exp_seq)
        mapped = set(self.mapped_id(seq_id) for seq_id in aligned)
        mapped.discard(None)

        if len(mapped) > 1:
            self.logger.error("Sequence %s has multiple possiblities %s %s:",
                              exp_seq, aligned, mapped)

        if len(mapped) != 1:
            return None
        return mapped.pop()

    def match(self, exp_seq, **kwargs):
        if 'filename' in kwargs:
            return self.ss_id(kwargs['filename'])
        return self.matching_sequence(exp_seq)

    def data(self, exp_seq, **kwargs):
        ss_id = self.match(exp_seq, **kwargs)
        if ss_id is None:
            raise core.Skip("Nothing to align %s to", exp_seq)

        return {
            'ss_id': ss_id,
            'exp_seq_id': exp_seq
        }
