import json

import core as core

from models import ExpSeqPosition
from models import ExpSeqInfo as ExpSeq
from models import CorrespondenceNts as Corr
from models import CorrespondenceInfo as Info
from models import PdbModifiedCorrespondecies

from utils.alignment import align
from utils.structures import NR as NrUtil
from utils.structures import Structure as StructureUtil


class Parser(object):
    def __init__(self, additional=None):
        self.additional = additional or {}

    def __call__(self, response):
        data = []
        for row in json.loads(response.text):
            entry = {'exp_seq_position_id1': row['unit1'],
                     'exp_seq_position_id2': row['unit2']}
            entry.update(self.additional)
            data.append(entry)
        return data


class Loader(core.Loader):
    name = 'correspondence_nts'
    update_gap = False

    def __init__(self, config, maker):
        self.nr_util = NrUtil(maker)
        self.st_util = StructureUtil(maker)
        super(Loader, self).__init__(config, maker)
        self.translation = self.__translation__()

    def structure_data(self, exp_id):
        ids = []
        sequence = []
        with self.session() as session:
            query = session.query(ExpSeqPosition).\
                filter(ExpSeqPosition.exp_seq_id == exp_id).\
                order_by(ExpSeqPosition.index)

            if not query.count():
                raise core.SkipValue("Could not get positions for %s" % exp_id)

            for result in query:
                seq_id = result.id
                seq = seq_id.split('|')[3]
                seq = self.translation.get(seq, seq)

                if seq not in ['A', 'C', 'G', 'U']:
                    raise core.SkipValue("Bad unit %s for %s" % seq, seq_id)

                ids.append(seq_id)
                sequence.append(seq)

        sequence = ''.join(sequence)
        return {'ids': ids, 'sequence': sequence}

    def correlate(self, corr_id, ref, target):
        if not ref or not ref.get('ids') or not ref.get('sequence'):
            self.logger.error("Not given complete reference: %s", ref)
            raise core.InvalidState("Incomplete reference")

        if not target or not target.get('ids') or not target.get('sequence'):
            self.logger.error("Not given complete target: %s", target)
            raise core.InvalidState("Incomplete target")

        results = align([ref, target])
        self.logger.info("Found %s correlations", len(results))
        data = []
        for result in results:
            data.append(Corr(exp_seq_position_id1=result[0],
                             exp_seq_position_id2=result[1],
                             correspondence_id=corr_id))
        return data

    def remove(self, corr_id):
        with self.session() as session:
            session.query(Corr).filter_by(correspondence_id=corr_id).\
                delete(synchronize_session=False)

    def has_data(self, corr_id):
        with self.session(log_exceptions=False) as session:
            query = session.query(Corr).filter_by(correspondence_id=corr_id)
            return bool(query.count())

    def transform(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(Info).\
                join(ExpSeq, ExpSeq.id == Info.exp_seq_id1).\
                filter(ExpSeq.pdb == pdb)
            return [result.id for result in query]

    def data(self, corr_id, **kwargs):
        with self.session() as session:
            query = session.query(Info).filter_by(id=corr_id)
            result = query.first()
            if not result:
                self.logger.error("Could not get corr with id %s", corr_id)
                raise core.InvalidState("Missing correspondence_id")
            exp1 = result.exp_seq_id1
            exp2 = result.exp_seq_id2

        self.logger.debug("Using sequences %s %s", exp1, exp2)
        exp1 = self.structure_data(exp1)
        exp2 = self.structure_data(exp2)

        return self.correlate(corr_id, exp1, exp2)

    def __translation__(self):
        mapping = {}
        with self.session() as session:
            query = session.query(PdbModifiedCorrespondecies)
            for result in query:
                mapping[result.modified_unit] = result.standard_unit
        return mapping
