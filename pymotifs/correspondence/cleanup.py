"""
This stage is present to clean up all alignments which are not good. We
store these alignments in the database temporarily, then summarize them. If
afterwards we see that these are not good alignments, we delete them from the
database to save space. Our selection criteria for alignments to do is very
broad so we end up aligning many things that aren't good matches to make sure
we do find all good matches. If we kept the all bad alignments we would waste
lots of space.
"""

from pymotifs import core
from pymotifs import models as mod

from pymotifs.correspondence.info import Loader as InfoLoader
from pymotifs.correspondence.summary import Loader as SummaryLoader


class Loader(core.Stage):
    mark = False
    dependencies = set([SummaryLoader, InfoLoader])

    def to_process(self, pdbs, **kwargs):
        """
        We transform the list of pdbs into the list of correspondences.

        :param list pdb: The list of pdb ids. Currently ignored.
        :param dict kwargs: The keyword arguments which are ignored.
        :returns: A list of correspondence ids to process.
        """

        if len(pdbs) < 500:
            raise core.Skip("Too few pdb files being processed to clean up correspondences")

        with self.session() as session:
            position = mod.CorrespondencePositions
            info = mod.CorrespondenceInfo
            subquery = session.query(position.correspondence_id).\
                distinct().\
                subquery()

            query = session.query(info).\
                join(subquery,
                     subquery.c.correspondence_id == info.correspondence_id)
            return [result.correspondence_id for result in query]

    def should_process(self, corr_id, **kwargs):
        """
        Test if we should process some correspondence. This is done by
        testing if it is not a good alignment. If not then we process,
        otherwise we do not.
        """

        with self.session() as session:
            corr = session.query(mod.CorrespondenceInfo).get(corr_id)
            return not bool(corr.good_alignment)

    def process(self, corr_id, **kwargs):
        """
        Process the correspodence. This will delete the correspondence positions
        as this means the correspondence is not a good alignment.
        """

        if kwargs.get('dry_run'):
            self.logger.info("Would positions for delete: %s", corr_id)
            return False

        with self.session() as session:
            positions = mod.CorrespondencePositions
            session.query(positions).\
                filter(positions.correspondence_id == corr_id).\
                delete(synchronize_session=False)
