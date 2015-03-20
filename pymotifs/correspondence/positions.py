from pymotifs import core as core

from pymotifs.models import ExpSeqInfo
from pymotifs.models import ExpSeqPosition
from pymotifs.models import ChainMapping
from pymotifs.models import ChainInfo
from pymotifs.models import CorrespondenceInfo as Info
from pymotifs.models import CorrespondencePositions as Position
from pymotifs.models import PdbModifiedCorrespondecies

from pymotifs.utils.alignment import align


class Loader(core.Loader):

    def __init__(self, config, maker):
        super(Loader, self).__init__(config, maker)
        self.translation = self.__translation__()
        self.valid_sequence = set(['A', 'C', 'G', 'U'])

    def remove(self, pdb):
        pass

    def has_data(self, pdb):
        pass

    def sequence(self, exp_id):
        ids = []
        sequence = []
        with self.session() as session:
            query = session.query(ExpSeqPosition).\
                filter(ExpSeqPosition.exp_seq_id == exp_id).\
                order_by(ExpSeqPosition.index)

            if not query.count():
                raise core.InvalidState("Could not get positions for %s" %
                                        exp_id)

            for result in query:
                seq_id = result.id
                seq = self.translation.get(result.unit, result.unit)

                if seq not in self.valid_sequence:
                    raise core.SkipValue("Bad unit %s for %s" % seq, seq_id)

                ids.append(seq_id)
                sequence.append(seq)

        sequence = ''.join(sequence)

        if not ids:
            raise core.InvalidState("Failed to get ids for %s" % exp_id)

        if not(sequence):
            raise core.InvalidState("Failed to get a sequence for %s" % exp_id)

        if len(ids) != len(sequence):
            raise core.InvalidState("Did not get as many ids as sequence")

        return {'ids': ids, 'sequence': sequence}

    def correlate(self, corr_id, ref, target):
        results = align([ref, target])
        self.logger.debug("Found %s correlations", len(results))
        data = []
        for result in results:
            data.append({
                'exp_seq_position_id1': result[0],
                'exp_seq_position_id2': result[1],
                'correspondence_id': corr_id
            })
        return data

    def __translation__(self):
        mapping = {}
        with self.session() as session:
            query = session.query(PdbModifiedCorrespondecies)
            for result in query:
                mapping[result.modified_unit] = result.standard_unit
        return mapping

    def pairs(self, pdb):
        pairs = []
        with self.session() as session:
            query = session.query(Info).\
                join(ExpSeqInfo,
                     ExpSeqInfo.id == Info.exp_seq_id1 |
                     ExpSeqInfo.id == Info.exp_seq_id2).\
                join(ChainMapping, ChainMapping.exp_seq_id == ExpSeqInfo.id).\
                join(ChainInfo, ChainInfo.id == ChainMapping.chain_id).\
                filter(ChainInfo.pdb_id == pdb)

            for result in query:
                pairs.append((result.id, result.exp_id1, result.exp_id2))

        return pairs

    def data(self, pdb, **kwargs):
        for corr_id, exp_id1, exp_id2 in self.pairs(pdb):
            sequence1 = self.sequence(exp_id1)
            sequence2 = self.sequence(exp_id2)
            for position in self.correlate(corr_id, sequence1, sequence2):
                yield Position(**position)
