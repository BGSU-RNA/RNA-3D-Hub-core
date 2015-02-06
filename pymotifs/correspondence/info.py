from pymotifs import core
from pymotifs.models import ExpSeqInfo as ExpSeq
from pymotifs.models import CorrespondenceInfo as Info

from pymotifs.utils.structures import NR as NrUtil
from pymotifs.utils.structures import Structure as StructureUtil


class MissingExpSeq(core.InvalidState):
    pass


class Loader(core.Loader):
    name = 'correspondence_info'
    update_gap = False
    allow_no_data = True
    stop_on_failure = False

    def __init__(self, config, maker):
        self.util = NrUtil(maker)
        self.structure = StructureUtil(maker)
        super(Loader, self).__init__(config, maker)

    def transform(self, pdb, **kwargs):
        """We transform the pdb into the exp_seq_id that corresponds to the
        longest chain in that structure. This means we will attempt to align
        only the longest chain to other similar structures.
        """
        exp_id = self.exp_seq_id(pdb)
        self.logger.info("Using exp_seq_id %s for %s", exp_id, pdb)
        return [exp_id]

    def exp_seq_id(self, pdb):
        chain = self.structure.longest_chain(pdb)

        with self.session() as session:
            query = session.query(ExpSeq).filter_by(pdb=pdb, chain=chain)
            result = query.first()
            if not result:
                self.logger.debug("No exp_seq for %s|%s", pdb, chain)
                return None
            result = result.id
        return result

    def pdb(self, exp_seq_id):
        with self.session() as session:
            query = session.query(ExpSeq).filter_by(id=exp_seq_id)
            result = query.first()
            if not result:
                self.logger.error("Could not find ExpSeq for given id")
                raise core.InvalidState("Missing exp_seq_id")
            return result.pdb, result.chain

    def possible(self, exp_seq_id):
        pdb, chain = self.pdb(exp_seq_id)
        members = set(self.util.members(pdb))

        possible = []
        for member in members:
            other_id = self.exp_seq_id(member)
            if other_id is None:
                self.logger.info("Skipping %s as it has no exp_seq_id", member)
            possible.append(other_id)
        return possible

    def known(self, exp_seq_id):
        with self.session() as session:
            query = session.query(Info).filter(Info.exp_seq_id1 == exp_seq_id)
            return [result.exp_seq_id2 for result in query]

    def missing(self, exp_seq_id):
        return set(self.possible(exp_seq_id)) - set(self.known(exp_seq_id))

    def has_data(self, exp_seq_id):
        """We check if we have all data for this exp_seq_id by seeing if we
        aligned this to all possible sequences.
        """
        return not bool(self.missing(exp_seq_id))

    def remove(self, pdb):
        """Removes all old correspondences for the given pdb. This will only
        remove the ones marked by possible.
        """
        possible = self.possible(pdb)
        with self.session() as session:
            session.query(Info).filter_by(pdb1=pdb).\
                filter(Info.pdb2.in_(possible)).\
                delete(synchronize_session=False)

    def data(self, exp_seq_id, **kwargs):
        data = []
        for other in self.missing(exp_seq_id):
            self.logger.debug("Using second exp_seq %s", other)
            data.append(Info(exp_seq_id1=exp_seq_id, exp_seq_id2=other))
        return data
