import sys
import logging
import traceback
import itertools as it

from MotifAtlasBaseClass import MotifAtlasBaseClass
import models as mod
import utils as ut

logger = logging.getLogger(__name__)


class MissingNucleotideException(Exception):
    """Raised when there is a missing nucleotide that we need to map.
    """
    pass


class StructureUtil(ut.DatabaseHelper):
    """Some useful methods
    """

    def loops(self, pdb):
        """Get all loops in a structure.
        """
        loops = []
        with self.session() as session:
            query = session.query(mod.LoopPositions).\
                join(mod.LoopsAll,
                     mod.LoopPositions.loop_id == mod.LoopsAll.id).\
                filter(mod.LoopsAll.pdb == pdb).\
                order_by(mod.LoopsAll.id)

            grouped = it.groupby(it.imap(ut.row2dict, query),
                                 lambda a: a['loop_id'])
            for loop_id, positions in grouped:
                loops.append({
                    'id': loop_id,
                    'nts': [pos['nt_id'] for pos in positions]
                })

        return loops

    def reference(self, pdb):
        """Get all correlated reference structures.
        """
        with self.session() as session:
            query = session.query(mod.PdbCorrespondences).filter_by(pdb2=pdb)
            return [result.pdb1 for result in query]

    def mapping(self, ref, pdb):
        """Get the mapping from nucleotides in the reference to the nucleotides
        in the pdb.
        """
    # PdbCorrespondences
        mapping = {}
        with self.session() as session:
            query = session.query(mod.NtNtCorrespondences).\
                join(mod.PdbCorrespondences,
                     mod.PdbCorrespondences.id ==
                     mod.NtNtCorrespondences.correspondence_id).\
                filter(mod.PdbCorrespondences.pdb1 == ref).\
                filter(mod.PdbCorrespondences.pdb2 == pdb)
            for result in query:
                mapping[result.unit1_id] = result.unit2_id
                mapping[result.unit2_id] = result.unit1_id
        return mapping


class Loader(MotifAtlasBaseClass, ut.DatabaseHelper):

    def __init__(self, maker):
        MotifAtlasBaseClass.__init__(self)
        ut.DatabaseHelper.__init__(self, maker)
        self._overlaps = {}

    def has_data(self, reference, pdb):
        with self.session() as session:
            query = self.__base_info_query__(session, reference, pdb)
            return bool(query.count())

    def remove_old(self, reference, pdb):
        pass

    def overlap(self, coverage):
        if not self._overlaps:
            with self.session() as session:
                for overlap in session.query(mod.Overlaps):
                    self._overlaps[overlap.name] = overlap.id

        if coverage not in self._overlaps:
            raise Exception("Unknown overlap name: " + coverage)

        return self._overlaps[coverage]

    def discrepancy(self, loop1, loop2):
        """Get the discrepancy between two loops.
        """
        return None

    def map(self, nts, mapping):
        for nt in nts:
            if nt not in mapping:
                raise MissingNucleotideException("Can't map missing nt: " + nt)
            yield mapping[nt]

    def coverage(self, loop1, loop2, mapping):
        """Get the coverage between the two loops.
        """
        loop_nt1 = set(loop1['nts'])

        try:
            mapped = set(self.map(loop2['nts'], mapping))
        except:
            logging.error("Could not map all nucleotides")
            logger.error(traceback.format_exc(sys.exc_info()))
            return None

        intersection = mapped.intersection(loop_nt1)
        if not intersection:
            return None
        if mapped == loop_nt1:
            return 'exact'
        if mapped.issubset(loop_nt1):
            return 'contained'
        if mapped.issuperset(loop_nt1):
            return 'enclose'
        if intersection:
            return 'partial'
        raise Exception("This should never occur")

    def compare(self, reference, pdb, mapping):
        loops = self.loops(pdb)

        for ref_loop in self.loops(reference):
            found = False
            for loop in loops:

                try:
                    cover = self.coverage(ref_loop, loop, mapping)
                except:
                    logging.error("Exception in overlap between %s, %s",
                                  ref_loop['id'], loop['id'])
                    logger.error(traceback.format_exc(sys.exc_info()))
                    continue

                disc = None
                if cover:
                    found = True
                    if cover == 'exact':
                        disc = self.discrepancy(ref_loop, loop)

                    yield mod.LoopOverlap(loop1_id=ref_loop['id'],
                                          loop2_id=loop['id'],
                                          overlap_id=self.overlap(cover),
                                          discrepancy=disc)

            if not found:
                    yield mod.LoopOverlap(loop1_id=ref_loop['id'],
                                          loop2_id=None,
                                          overlap_id=self.overlap('unique'),
                                          discrepancy=None)

    def data(self, pdb, recalculate=False, **kwargs):
        for ref_pdb in self.reference(pdb):
            mapping = self.mapping(ref_pdb, pdb)
            yield it.chain(self.compare(pdb, ref_pdb, mapping),
                           self.compare(ref_pdb, pdb, mapping))

    def __call__(self, pdbs, **kwargs):
        if not pdbs:
            raise Exception("No pdbs given")

        for pdb in pdbs:
            logger.info("Getting loop loop correspondence for %s", pdb)

            try:
                data = self.data(pdb, **kwargs)
                self.store(data)
            except:
                logger.error("Failed storing loop loop correspondencies in %s",
                             pdb)
                logger.error(traceback.format_exc(sys.exc_info()))


if __name__ == '__main__':
    from utils import main
    main(Loader)
