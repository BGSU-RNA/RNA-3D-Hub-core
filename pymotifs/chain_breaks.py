import logging

from models import PolymerInfo
import utils as ut
import core

from rnastructure.tertiary.cif import CIF
from rnastructure.util import unit_ids as uids

logger = logging.getLogger(__name__)


class ChainBreakFinder(object):

    def __call__(self, cif_file):
        with open(cif_file, 'rb') as raw:
            cif = CIF(raw)

        breaks = []
        for polymer in cif.polymers():
            breaks.append((polymer.first().unit_id(),
                           polymer.last().unit_id()))
        return breaks


class Loader(core.Loader):
    finder = ChainBreakFinder()

    def __init__(self, config, maker):
        super(Loader, self).__init__(self, config, maker)
        self.cif = ut.CifFileFinder(self.config)

    def has_data(self, pdb):
        with self.session() as session:
            query = session.query(PolymerInfo).\
                filter_by(pdb_id=pdb)
            return bool(query.count())

    def remove(self, pdb):
        with self.session() as session:
            session.query(PolymerInfo).filter_by(pdb_id=pdb).delete()

    def data(self, pdb, **kwargs):
        cif_file = self.cif(pdb)
        endpoints = self.finder(cif_file)
        converter = uids.generate_converter('unit', 'nucleotide')

        data = []
        for unit1_id, unit2_id in endpoints:
            parts = unit1_id.split('|')
            # TODO: Do not assume AU only
            nt1_id = converter(unit1_id, type='AU')
            nt2_id = converter(unit2_id, type='AU')
            data.append(PolymerInfo(start_unit_id=nt1_id,
                                    end_unit_id=nt2_id,
                                    chain=parts[2],
                                    model=int(parts[1]),
                                    pdb_id=pdb))
        return data


if __name__ == '__main__':
    from utils import main
    main(Loader)
