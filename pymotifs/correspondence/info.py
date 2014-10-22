import logging

import core
from models import CorrespondenceInfo as Info

from utils.structures import NR as NrUtil

logger = logging.getLogger(__name__)


class Loader(core.Loader):
    name = 'correspondence_info'
    update_gap = False

    def __init__(self, config, maker):
        self.util = NrUtil(maker)
        super(Loader, self).__init__(config, maker)

    def possible(self, pdb):
        """Get all possible pdbs to align against. This means all pdbs in the
        same nr class as the given pdb for the current release.
        """
        return self.util.members(pdb)

    def known(self, pdb):
        """Get all pdbs we have aligned against.
        """
        with self.session() as session:
            query = session.query(Info).filter(Info.pdb1 == pdb)
            return [result.pdb2 for result in query]

    def missing(self, pdb):
        """Determine which pdbs we have not yet aligned against.
        """
        return set(self.possible(pdb)) - set(self.known(pdb))

    def has_data(self, pdb):
        """We check if we have all data for this pdb by seeing if we aligned
        this to all possible pdbs.
        """
        return not bool(self.missing(pdb))

    def remove(self, pdb):
        """Removes all old correspondences for the given pdb. This will only
        remove the ones marked by possible.
        """
        logging.info("Clearing out all old correspondence_info for %s", pdb)
        possible = self.possible(pdb)
        with self.session() as session:
            session.query(Info).filter_by(pdb1=pdb).\
                filter(Info.pdb2.in_(possible)).\
                delete(synchronize_session=False)

    def data(self, pdb, **kwargs):
        return [Info(pdb1=pdb, pdb2=other) for other in self.possible(pdb)]


if __name__ == '__main__':
    from utils import main
    main(Loader)
