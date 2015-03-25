"""This stage is present to clean up all alignments which are not good. We
store these alignments in the database temporarily, then summarize them. If
afterwards we see that these are not good alignments, we delete them from the
database to save space. Our selection criteria for alignments to do is very
broad so we end up aligning many things that aren't good matches to make sure
we do find all good matches. If we kept the all bad alignments we would waste
lots of space.
"""

from sqlalchemy.orm import aliased

from pymotifs import core
from pymotifs import utils

from pymotifs.models import ExpSeqInfo as ExpInfo
from pymotifs.models import CorrespondenceInfo as Info
from pymotifs.models import CorrespondencePositions as Position


class Loader(core.Stage):
    mark = False

    def info(self, corr_id):
        with self.session() as session:
            info = session.query(Info).get(corr_id)
            info = utils.row2dict(info)

        with self.session() as session:
            e1 = aliased(ExpInfo)
            e2 = aliased(ExpInfo)

            query = session.query(Info.id, e1.length.label('first'),
                                  e2.length.label('second')).\
                join(e1, Info.exp_seq_id1 == e1.id).\
                join(e2, Info.exp_seq_id2 == e2.id).\
                filter(Info.id == corr_id)
            result = query.one()

            info['min'] = min(result.first, result.second)

        return info

    def should_process(self, corr_id, **kwargs):
        """Detect if the given correspondence id is below our cutoffs for a
        good match.
        """

        info = self.info(corr_id)
        if info['min'] <= 36:
            return not bool(info['mismatch_count'])
        if info['min'] <= 80:
            return info['mismatch_count'] > 4
        return 0.9 <= (float(info['mismatch_count']) / float(info['min']))

    def to_process(self, pdbs, **kwargs):
        """We transform the list of pdbs into the list of correspondences.

        :param list pdb: The list of pdb ids. Currently ignored.
        :param dict kwargs: The keyword arguments which are ignored.
        :returns: A list of correspondence ids to process.
        """

        with self.session() as session:
            query = session.query(Info)
            return [result.id for result in query]

    def process(self, corr_id, **kwargs):
        """Process the correspodence. This will delete the correspondence positions
        as this is not a good alignment.
        """

        if kwargs.get('dry_run'):
            self.logger.info("Would positions for delete: %s", corr_id)
            return False

        with self.session() as session:
            session.query(Position).\
                filter(Position.correspondence_id == corr_id).\
                delete(synchronize_session=False)
