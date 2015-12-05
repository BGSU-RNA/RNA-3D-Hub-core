"""This is a stage to align two experimental sequences and store the alignment
between each nucleotide.
"""

from pymotifs import core as core

from pymotifs.models import ExpSeqPosition
from pymotifs.models import CorrespondenceInfo as Info
from pymotifs.models import CorrespondencePositions as Position

from pymotifs.correspondence.info import Loader as CorrLoader
from pymotifs.exp_seq.info import Loader as InfoLoader
from pymotifs.exp_seq.positions import Loader as PositionLoader

from pymotifs.utils.alignment import align


class Loader(core.SimpleLoader):
    """A loader for computing the position to position alignment and storing
    it.
    """

    mark = False
    dependencies = set([CorrLoader, InfoLoader, PositionLoader])

    def to_process(self, pdbs, **kwargs):
        """We transform all the pdbs into the correspodencies to do. While this
        does not respect the pdbs given but it does make all other code a lot
        cleaner and easier to understand.

        :param list pdb: The list of pdb ids. Currently ignored.
        :param dict kwargs: The keyword arguments which are ignored.
        :returns: A list of correspodence ids to process.
        """

        with self.session() as session:
            query = session.query(Info).filter(Info.length != None)
            return [result.correspondence_id for result in query]

    def query(self, session, corr_id):
        with self.session() as session:
            return session.query(Position).filter_by(correspondence_id=corr_id)

    def sequence(self, exp_id):
        """Load all information about the experimental sequence with the given
        id. This will load both the ids and the sequence. The ids will be a
        list of numbers, while the sequence is a string.
        """

        ids = []
        sequence = []
        with self.session() as session:
            query = session.query(ExpSeqPosition).\
                filter(ExpSeqPosition.exp_seq_id == exp_id).\
                order_by(ExpSeqPosition.index)

            if not query.count():
                self.logger.error("Could not get sequence for %s", exp_id)
                raise core.Skip("Skipping correspondence")

            for index, result in enumerate(query):
                seq_id = result.exp_seq_position_id
                seq = result.normalized_unit or 'N'
                ids.append(seq_id)
                sequence.append(seq)

        return {'ids': ids, 'sequence': ''.join(sequence)}

    def correlate(self, corr_id, ref, target):
        results = align([ref, target])
        self.logger.debug("Alignment is %i long", len(results))
        data = []
        for index, result in enumerate(results):
            data.append({
                'exp_seq_position_id_1': result[0],
                'exp_seq_position_id_2': result[1],
                'correspondence_id': corr_id,
                'index': index
            })
        return data

    def info(self, corr_id):
        with self.session() as session:
            query = session.query(Info).filter_by(correspondence_id=corr_id)
            result = query.one()
            return (result.exp_seq_id_1, result.exp_seq_id_2)

    def data(self, corr_id, **kwargs):
        exp_id1, exp_id2 = self.info(corr_id)
        sequence1 = self.sequence(exp_id1)
        sequence2 = self.sequence(exp_id2)

        for position in self.correlate(corr_id, sequence1, sequence2):
            yield Position(**position)

            pos1 = position['exp_seq_position_id_1']
            pos2 = position['exp_seq_position_id_2']
            if pos1 != pos2:
                rev = dict(position)
                rev['exp_seq_position_id_1'] = pos2
                rev['exp_seq_position_id_2'] = pos1

                yield Position(**rev)
