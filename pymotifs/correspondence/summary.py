"""This is a loader to summarize the correspondences. This requires that the
positions have already been loaded, or bad things will happen. The
summarization is key to building the nr sets correctly.
"""

from sqlalchemy.orm import aliased

from pymotifs import core
from pymotifs import utils
from pymotifs import models as mod

from pymotifs.correspondence.positions import Loader as PositionLoader

from pymotifs.constants import CORRESPONDENCE_EXACT_CUTOFF
from pymotifs.constants import CORRESPONDENCE_LIMITED_CHANGES



class Loader(core.Loader):
    merge_data = True
    mark = False

    dependencies = set([PositionLoader])

    def to_process(self, pdbs, **kwargs):
        """We transform the list of pdbs into the list of correspondences that
        have not yet been summarized.

        :param list pdb: The list of pdb ids. Currently ignored.
        :param dict kwargs: The keyword arguments which are ignored.
        :returns: A list of correspondence ids to process.
        """

        with self.session() as session:
            query = session.query(mod.CorrespondenceInfo.correspondence_id)
            return [result.correspondence_id for result in query]

    def remove(self, corr_id, **kwargs):
        """We do not remove anything when summarizing as we aren't actually
        adding any rows when we do this, to force a recompute of this data you
        should do a recompute on the whole correspondence level information.
        """
        self.logger.info("Not removing anything, recompute all correspondence")

    def has_data(self, corr_id, **kwargs):
        """Check if we have summarized this correspondence before. This only
        looks for the length field not being null.
        """

        with self.session() as session:
            query = session.query(mod.CorrespondenceInfo).\
                filter(mod.CorrespondenceInfo.correspondence_id == corr_id).\
                filter(mod.CorrespondenceInfo.length != None)
            return bool(query.count())

    def current(self, corr_id):
        """Get the current data for the correspondence.
        """

        with self.session() as session:
            info = session.query(mod.CorrespondenceInfo).get(corr_id)
            return utils.row2dict(info)                                     ## A dict for the values of columns in correspondence_info table for a specific corr_id.

    def sizes(self, info):
        """Compute the minimum size of the experimental sequences used in this
        correspondence.

        :param dict info: The information about the correspondence
        :returns: The minimum size
        """

        with self.session() as session:
            e1 = aliased(mod.ExpSeqInfo)
            e2 = aliased(mod.ExpSeqInfo)
            corr = info['correspondence_id']

            query = session.query(mod.CorrespondenceInfo.correspondence_id,
                                  e1.length.label('first'),
                                  e2.length.label('second')).\
                join(e1,
                     mod.CorrespondenceInfo.exp_seq_id_1 == e1.exp_seq_id).\
                join(e2,
                     mod.CorrespondenceInfo.exp_seq_id_2 == e2.exp_seq_id).\
                filter(mod.CorrespondenceInfo.correspondence_id == corr)
            result = query.one()

            return sorted([result.first, result.second])

    def good_alignment(self, info, min_size, max_size, **kwargs):
        """Detect if the given correspondence id is below our cutoffs for a
        good match.
        """

        if not info['aligned_count']:
            return False

        if min_size < CORRESPONDENCE_EXACT_CUTOFF:                           ## ====> min_size < 19
            if min_size == max_size:
                return info['mismatch_count'] == 0
            return False

        if min_size < CORRESPONDENCE_LIMITED_CHANGES:                        ## ===> if min_size < 80
            return info['mismatch_count'] <= 4

        if max_size > min_size * 2:               
            return False

        if not info['mismatch_count']:
            return True

        return float(info['match_count']) / float(min_size) >= 0.95          

    def alignment(self, corr_id):
        with self.session() as session:
            p1 = aliased(mod.ExpSeqPosition)
            p2 = aliased(mod.ExpSeqPosition)
            query = session.query(mod.CorrespondencePositions.correspondence_positions_id,
                                  p1.unit.label('unit1'),
                                  p2.unit.label('unit2')).\
                filter(mod.CorrespondencePositions.correspondence_id == corr_id).\
                outerjoin(p1,
                          p1.exp_seq_position_id == mod.CorrespondencePositions.exp_seq_position_id_1).\
                outerjoin(p2,
                          p2.exp_seq_position_id == mod.CorrespondencePositions.exp_seq_position_id_2).\
                order_by(mod.CorrespondencePositions.index).\
                group_by(mod.CorrespondencePositions.index)

            results = []
            for result in query:
                results.append({'unit1': result.unit1, 'unit2': result.unit2}) 
        return results

    def summary(self, positions):
        data = {
            'length': len(positions),
            'aligned_count': 0,
            'first_gap_count': 0,
            'second_gap_count': 0,
            'match_count': 0,
            'mismatch_count': 0
        }

        seq1 = ''
        seq2 = ''
        method2_seq1 = []
        method2_seq2 = []

        for position in positions:                      ## the following part is for matching bases for two sequences
            unit1 = position['unit1']
            unit2 = position['unit2']

            seq1 = seq1+unit1
            seq2 = seq2+unit2
            method2_seq1.append(unit1)
            method2_seq2.append(unit2)

            if unit1 and unit2:
                data['aligned_count'] += 1
                if unit1 == unit2:
                    data['match_count'] += 1
                else:
                    data['mismatch_count'] += 1
            else:
                data['mismatch_count'] += 1

            if not unit1:
                data['first_gap_count'] += 1
            if not unit2:
                data['second_gap_count'] += 1

        return data

    def data(self, corr_id, **kwargs):
        """Compute the summary for the given correspondence id. This will
        update the entry with the counts of match, mismatch and such.
        """

        data = self.current(corr_id)
        min_size, max_size = self.sizes(data)
        data.update(self.summary(self.alignment(corr_id)))
        data['good_alignment'] = self.good_alignment(data, min_size, max_size)

        return mod.CorrespondenceInfo(**data)
