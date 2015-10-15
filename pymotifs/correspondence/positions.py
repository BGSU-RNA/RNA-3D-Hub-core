"""This is a stage to align two experimental sequences and store the alignment
between each nucleotide.
"""

from pymotifs import core as core

from pymotifs.models import ExpSeqPosition
from pymotifs.models import CorrespondenceInfo as Info
from pymotifs.models import CorrespondencePositions as Position
from pymotifs.models import RnaUnitModifiedCorrespondencies
from pymotifs.exp_seq.info import Loader as InfoLoader
from pymotifs.exp_seq.positions import Loader as PositionLoader

from pymotifs.utils.alignment import align


class Loader(core.Loader):
    """A loader for computing the position to position alignment and storing
    it.
    """
    mark = False

    dependencies = set([InfoLoader, PositionLoader])

    def __init__(self, config, maker):
        super(Loader, self).__init__(config, maker)
        self.translation = self.__translation__()
        self.valid_sequence = set(['A', 'C', 'G', 'U', 'N'])

    def to_process(self, pdbs, **kwargs):
        """We transform all the pdbs into the correspodencies to do. While this
        does not respect the pdbs given but it does make all other code a lot
        cleaner and easier to understand.

        :param list pdb: The list of pdb ids. Currently ignored.
        :param dict kwargs: The keyword arguments which are ignored.
        :returns: A list of correspodence ids to process.
        """

        with self.session() as session:
            query = session.query(Info)
            return [result.correspondence_id for result in query]

    def has_data(self, corr_id, **kwargs):
        """This is an unusual check for data. We do not first look up to see if
        positions with the given corr_id exist. If that is not the case, we
        then check to see if in the info table we have summerized the
        alignment. In that case we do not need to recompute the alignment and
        store the positions. We later cleanup the alignments to remove any bad
        ones, so we cannot just check that we need to align by using the
        positions.
        """

        with self.session() as session:
            query = session.query(Position).\
                filter(Position.correspondence_id == corr_id)

            if query.count():
                return True

        with self.session() as session:
            query = session.query(Info).\
                filter(Info.length != None).\
                filter(Info.correspondence_id == corr_id)
            return bool(query.count())

    def remove(self, corr_id, **kwargs):
        with self.session() as session:
            session.query(Position).\
                filter(Position.correspondence_id == corr_id).\
                delete(synchronize_session=False)

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
                return None

            size = int(query.count()) - 1
            for index, result in enumerate(query):
                seq_id = result.id
                seq = self.translation.get(result.unit, result.unit)

                # X is really N, but old PDB data doesn't respect that.
                if 'X' in seq:
                    self.logger.debug("For sequence %s, changing X to N",
                                      exp_id)
                    seq = seq.replace('X', 'N')

                if seq not in self.valid_sequence:
                    if index != size:
                        raise core.Skip("Bad unit %s for %s" % (seq, seq_id))

                    self.logger.warning(("Skipping bad unit %s for %s as"
                                         "it may be a tRNA/AA"), seq, seq_id)
                    break

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
            return (result.exp_seq_id1, result.exp_seq_id2)

    def data(self, corr_id, **kwargs):
        exp_id1, exp_id2 = self.info(corr_id)
        sequence1 = self.sequence(exp_id1)
        sequence2 = self.sequence(exp_id2)

        if not sequence1:
            self.logger.error(("Could not get positions for experimental "
                               "sequence %s") % exp_id1)

        if not sequence2:
            self.logger.error(("Could not get positions for experimental "
                               "sequence %s") % exp_id2)

        if not sequence1 or not sequence2:
            raise core.Skip("Skipping this correspodence %s" % corr_id)

        for position in self.correlate(corr_id, sequence1, sequence2):
            yield Position(**position)

            pos1 = position['exp_seq_position_id_1']
            pos2 = position['exp_seq_position_id_2']
            if pos1 != pos2:
                rev = dict(position)
                rev['exp_seq_position_id_1'] = pos2
                rev['exp_seq_position_id_2'] = pos1

                yield Position(**rev)

    def __translation__(self):
        mapping = {}
        with self.session() as session:
            query = session.query(RnaUnitModifiedCorrespondencies)
            for result in query:
                mapping[result.id] = result.standard_unit
        return mapping
