import sys
import logging

from MotifAtlasBaseClass import MotifAtlasBaseClass
from models import PolymerInfo
from utils import DatabaseHelper
from utils import CifFileFinder

sys.path.append('rnastructure')
from rnastructure.tertiary.cif import CIF
from rnastructure.util import unit_ids as uids

# logger = logging.getLogger('chain_breaks')
logger = logging


class ChainBreakFinder(object):

    def __call__(self, cif_file):
        with open(cif_file, 'rb') as raw:
            cif = CIF(raw)

        breaks = []
        for polymer in cif.polymers():
            breaks.append((polymer.first().unit_id(),
                           polymer.last().unit_id()))
        return breaks


class ChainBreakLoader(MotifAtlasBaseClass, DatabaseHelper):
    finder = ChainBreakFinder()

    def __init__(self, maker):
        MotifAtlasBaseClass.__init__(self)
        self.cif = CifFileFinder(self.config)
        DatabaseHelper.__init__(self, maker)

    def has_breaks(self, pdb):
        with self.session() as session:
            query = session.query(PolymerInfo).\
                filter_by(pdb_id=pdb)
            return bool(query.count())

    def data(self, pdb):
        cif_file = self.cif(pdb)
        endpoints = self.finder(cif_file)
        converter = uids.generate_converter('unit', 'nucleotide')

        data = []
        for unit1_id, unit2_id in endpoints:
            parts = unit1_id.split('|')
            nt1_id = converter(unit1_id, type='AU')
            nt2_id = converter(unit2_id, type='AU')
            data.append(PolymerInfo(start_unit_id=nt1_id,
                                    end_unit_id=nt2_id,
                                    chain=parts[2],
                                    model=int(parts[1]),
                                    pdb_id=pdb))
        return data

    def __call__(self, pdbs):
        if not pdbs:
            raise Exception("No pdbs given")

        for pdb in pdbs:
            if self.has_breaks(pdb):
                logger.info("Skipping getting breaks for %s", pdb)
                continue

            logger.info("Getting breaks for %s", pdb)
            breaks = self.data(pdb)
            logger.info("Found %s breaks", len(breaks))

            try:
                logger.debug("Storing breaks for %s", pdb)
                self.store(breaks)
            except:
                logger.error("Failed to store breaks for %s", pdb)


if __name__ == '__main__':
    from utils import main
    main(ChainBreakLoader)
