import itertools as it
import collections as coll

import core

from models import ExpSeqInfo as ExpSeq
from models import CorrespondenceInteractions as CorrInts
from models import CorrespondenceInfo as Info

from utils.structures import Structure

QUERY = '''
select *
from correspondence_units
where
    correspondence_id = :val
    and old_id2 not like '%_BA1_%'
;
'''


class Loader(core.Loader):
    name = 'correspondence_interactions'
    update_gap = False

    def __init__(self, config, maker):
        self.util = Structure(maker)
        super(Loader, self).__init__(config, maker)

    def load_exp_info(self, corr_id, exp_seq_id):
        column = getattr(Info, exp_seq_id)
        with self.session() as session:
            query = session.query(ExpSeq).\
                join(Info, column == ExpSeq.id).\
                filter(Info.id == corr_id)

            result = query.first()
            if not result:
                raise core.SkipValue("Cannot get exp_seq_id: %s" % exp_seq_id)

            return result.pdb, result.chain

    def load_interaction(self, corr_id, exp_seq):
        pdb, chain = self.load_exp_info(corr_id, exp_seq)
        return self.util.interactions(pdb, chain)

    def nt_pairs(self, interactions):
        pairs = it.imap(lambda inter: (inter['iPdbSig'], inter['jPdbSig']),
                        interactions)
        return pairs

    def pair_mapping(self, interactions):
        pairs = self.nt_pairs(interactions)
        ids = it.imap(lambda inter: inter['id'], interactions)
        zipped = it.izip(pairs, ids)
        return dict(zipped)

    def translation(self, corr_id):
        mapping = {}
        with self.session() as session:
            result = session.execute(QUERY, {'val': corr_id})
            Record = coll.namedtuple('Record', result.keys())
            for result in it.imap(lambda r: Record(*r), result.fetchall()):
                mapping[result.old_id1] = result.old_id2
        return mapping

    def transform(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(Info).\
                join(ExpSeq, ExpSeq.id == Info.exp_seq_id1).\
                filter(ExpSeq.pdb == pdb)
            return [result.id for result in query]

    def has_data(self, corr_id):
        with self.session() as session:
            query = session.query(CorrInts).\
                filter_by(correspondence_id=corr_id)
            return bool(query.count())

    def remove(self, corr_id):
        with self.session() as session:
            session.query(CorrInts).\
                filter_by(correspondence_id=corr_id).\
                delete(synchronize_session=False)

    def data(self, corr_id, **kwargs):
        interactions = self.load_interaction(corr_id, 'exp_seq_id1')
        second = self.load_interaction(corr_id, 'exp_seq_id2')
        pair_mapping = self.pair_mapping(second)
        translation = self.translation(corr_id)
        pairs = self.nt_pairs(interactions)
        ids = it.imap(lambda inter: inter['id'], interactions)

        data = []
        for id1, (nt1, nt2) in it.izip(ids, pairs):
            trans1 = translation.get(nt1)
            trans2 = translation.get(nt2)
            if not trans1 or not trans2:
                continue

            trans = (trans1, trans2)
            id2 = pair_mapping.pop(trans, None)

            data.append(CorrInts(interaction_id1=id1, interaction_id2=id2,
                                 correspondence_id=corr_id))

        for id in pair_mapping.values():
            data.append(CorrInts(interaction_id1=None, interaction_id2=id,
                                 correspondence_id=corr_id))

        return data
