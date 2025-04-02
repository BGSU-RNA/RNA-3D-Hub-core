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

        # if len(pdbs) < 500:
        #     raise core.Skip("Too few pdb files being processed to clean up correspondences")

        with self.session() as session:
            # get correspondence_id values that are in both tables
            position = mod.CorrespondencePositions
            subquery = session.query(position.correspondence_id).\
                distinct().\
                subquery()

            info = mod.CorrespondenceInfo
            query = session.query(info).\
                join(subquery,
                     subquery.c.correspondence_id == info.correspondence_id).\
                filter(info.good_alignment == False)

            poor_alignment_ids = set([result.correspondence_id for result in query])

        if len(poor_alignment_ids) == 0:
            raise core.Skip("No poor alignments found")

        return poor_alignment_ids

    def should_process(self, corr_id, **kwargs):
        """
        Test if we should process some correspondence. This is done by
        testing if it is not a good alignment. If not then we process,
        otherwise we do not.

        It's slow to check this separately for each correspondence_id
        identified in to_process, but faster now that to_process only
        returns those correspondence_id values that have a poor alignment.
        Once we are confident that it's doing OK, we could have this
        function always return True, that would be faster.
        """

        with self.session() as session:
            corr = session.query(mod.CorrespondenceInfo).get(corr_id)
            return not bool(corr.good_alignment)

    def process(self, corr_id, **kwargs):
        """
        Process the correspondence.
        This will delete the correspondence positions
        because the correspondence is not a good alignment.
        But it will leave the entry in correspondence.info so we can
        see that we tried.
        """

        if kwargs.get('dry_run'):
            self.logger.info("Would positions for delete: %s", corr_id)
            return False

        with self.session() as session:
            positions = mod.CorrespondencePositions
            session.query(positions).\
                filter(positions.correspondence_id == corr_id).\
                delete(synchronize_session=False)
