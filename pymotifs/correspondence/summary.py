"""This is a loader to summarize the correspondences. This requires that the
positions have already been loaded, or bad things will happen. The
summarization is key to building the nr sets correctly.
"""

from sqlalchemy.orm import aliased

from pymotifs import core
from pymotifs import utils

from pymotifs.models import ExpSeqPosition as ExpPosition
from pymotifs.models import CorrespondenceInfo as Info
from pymotifs.models import CorrespondencePositions as Position


class Loader(core.Loader):
    merge_data = True
    mark = False

    def to_process(self, pdbs, **kwargs):
        """We transform the list of pdbs into the list of correspondences that
        have not yet been summarized.

        :param list pdb: The list of pdb ids. Currently ignored.
        :param dict kwargs: The keyword arguments which are ignored.
        :returns: A list of correspondence ids to process.
        """

        with self.session() as session:
            query = session.query(Info)
            return [result.id for result in query]

    def remove(self, corr_id, **kwargs):
        """We do not remove anything when summarizing as we aren't actually
        adding any rows when we do this, to force a recompute of this data you
        should do a recompute on the whole correspondence level information.
        """
        pass

    def has_data(self, corr_id, **kwargs):
        """Check if we have summarized this correspondence before. This only
        looks for the length field not being null.
        """

        with self.session() as session:
            query = session.query(Info).\
                filter(Info.id == corr_id).\
                filter(Info.length != None)
            return bool(query.count())

    def current(self, corr_id):
        """Get the current data for the correspondence.
        """

        with self.session() as session:
            info = session.query(Info).get(corr_id)
            return utils.row2dict(info)

    def data(self, corr_id, **kwargs):
        """Compute the summary for the given correspondence id. This will
        update the entry with the counts of match, mismatch and such.
        """

        data = self.current(corr_id)
        data.update({
            'length': 0,
            'aligned_count': 0,
            'first_gap_count': 0,
            'second_gap_count': 0,
            'match_count': 0,
            'mismatch_count': 0
        })

        with self.session() as session:
            p1 = aliased(ExpPosition)
            p2 = aliased(ExpPosition)
            query = session.query(Position.id,
                                  p1.unit.label('unit1'),
                                  p2.unit.label('unit2')).\
                filter(Position.correspondence_id == corr_id).\
                outerjoin(p1, p1.id == Position.exp_seq_position_id1).\
                outerjoin(p2, p2.id == Position.exp_seq_position_id2).\
                order_by(Position.index)

            data['length'] = query.count()
            for result in query:
                unit1 = result.unit1
                unit2 = result.unit2

                if unit1 and unit2:
                    data['aligned_count'] += 1
                    if unit1 == unit2:
                        data['match_count'] += 1
                    else:
                        data['mismatch_count'] += 1
                if not unit1:
                    data['first_gap_count'] += 1
                if not unit2:
                    data['second_gap_count'] += 1

        return Info(**data)
