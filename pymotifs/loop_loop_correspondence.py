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


class InvalidCoverageState(Exception):
    """Raised when something happens to put is in an invalid state when getting
    the coverage.
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
            query = session.query(mod.CorrespondenceInfo).filter_by(pdb2=pdb)
            return [ut.row2dict(result) for result in query]

    def mapping(self, ref, pdb):
        """Get the mapping from nucleotides in the reference to the nucleotides
        in the pdb.
        """
        mapping = {}
        with self.session() as session:
            query = session.query(mod.CorrespondenceNts).\
                join(mod.CorrespondenceInfo,
                     mod.CorrespondenceInfo.id ==
                     mod.CorrespondenceNts.correspondence_id).\
                filter(mod.CorrespondenceInfo.pdb1 == ref).\
                filter(mod.CorrespondenceInfo.pdb2 == pdb)
            for result in query:
                mapping[result.unit1_id] = result.unit2_id
                mapping[result.unit2_id] = result.unit1_id

        if not mapping:
            logger.error("Could not generate mapping between %s to %s", ref,
                         pdb)

        return mapping


class Loader(MotifAtlasBaseClass, ut.DatabaseHelper):

    def __init__(self, maker):
        MotifAtlasBaseClass.__init__(self)
        ut.DatabaseHelper.__init__(self, maker)
        self._overlaps = {}
        self.utils = StructureUtil(maker)

    def has_data(self, reference, pdb):
        with self.session() as session:
            query = self.__base_info_query__(session, reference, pdb)
            return bool(query.count())

    def remove_old(self, reference, pdb):
        pass

    def overlap(self, coverage):
        if not self._overlaps:
            with self.session() as session:
                for overlap in session.query(mod.LoopOverlapInfo):
                    self._overlaps[overlap.name] = overlap.id

        if coverage not in self._overlaps:
            raise Exception("Unknown overlap name: " + coverage)

        return self._overlaps[coverage]

    def discrepancy(self, loop1, loop2):
        """Get the discrepancy between two loops.
        """
        return None

    def map(self, nts, mapping):
        if not mapping:
            logger.error("Given empty mapping. Invalid state")
            raise MissingNucleotideException("Empty mapping")

        for nt in nts:
            if nt not in mapping:
                raise MissingNucleotideException(nt)
            yield mapping[nt]

    def map_loops(self, loops, ref, mapping):
        mapped = []
        for loop in loops:
            try:
                nts = list(self.map(loop['nts'], mapping))
            except MissingNucleotideException as err:
                logger.warn("Removing loop %s which cannot be mapped to %s",
                            loop['id'], ref)
                logger.warn("Missing nt %s", str(err))
                continue

            mapped.append({'id': loop['id'], 'nts': nts})

        return mapped

    def coverage(self, loop1, loop2):
        """Get the coverage between the two loops.
        """
        if not loop1:
            raise InvalidCoverageState("Empty first loop")

        if not loop2:
            raise InvalidCoverageState("Empty second loop")

        loop_nt1 = set(loop1.get('nts', []))
        if not loop_nt1:
            raise InvalidCoverageState("First loop had no nts")

        nts = loop2.get('nts', [])
        if not nts:
            raise InvalidCoverageState("Second loop has no nts")

        mapped = set(loop2['nts'])
        if len(mapped) != len(loop2['nts']):
            raise InvalidCoverageState("Could not map all loop2 nts")

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

    def loop_comparision_id(self, loop1, loop2, discrepancy):
        with self.session() as session:
            query = session.query(mod.LoopLoopComparisions).\
                filter_by(loop1_id=loop1, loop2_id=loop2)

            if query.count() == 0:
                session.add(mod.LoopLoopComparisions(loop1_id=loop1,
                                                     loop2_id=loop2,
                                                     discrepancy=discrepancy))
                session.commit()

        with self.session() as session:
            return session.query(mod.LoopLoopComparisions).\
                filter_by(loop1_id=loop1, loop2_id=loop2).\
                one().id

    def compare_loop(self, ref_loop, loops, corr_id):
        overlapping = []
        for loop in loops:
            try:
                cover = self.coverage(ref_loop, loop)
            except:
                logger.error("Exception in overlap between %s, %s",
                             ref_loop['id'], loop['id'])
                logger.error(traceback.format_exc(sys.exc_info()))
                continue

            if not cover:
                continue

            disc = None
            if cover == 'exact':
                disc = self.discrepancy(ref_loop, loop)

            compare_id = self.loop_comparision_id(ref_loop['id'], loop['id'],
                                                  disc)

            overlapping.append((loop['id'],
                               mod.CorrespondenceLoops(
                                   loop_loop_comparisions_id=compare_id,
                                   correspondence_id=corr_id,
                                   loop_overlap_info_id=self.overlap(cover))))

        return overlapping

    def compare(self, ref_loops, loops, corr_id):
        unseen_loops = set(loop['id'] for loop in loops)
        for ref in ref_loops:
            overlapping = self.compare_loop(ref, loops, corr_id)
            unseen_loops -= set(loop[0] for loop in overlapping)

            if overlapping:
                for compare in overlapping:
                    yield compare[1]

            else:
                compare_id = self.loop_comparision_id(ref['id'], None, None)
                yield mod.CorrespondenceLoops(
                    loop_loop_comparisions_id=compare_id,
                    correspondence_id=corr_id,
                    loop_overlap_info_id=self.overlap('unique'))

        for loop in unseen_loops:
            compare_id = self.loop_comparision_id(None, loop, None)

            yield mod.CorrespondenceLoops(
                loop_loop_comparisions_id=compare_id,
                correspondence_id=corr_id,
                loop_overlap_info_id=self.overlap('unique'))

    def data(self, pdb, recalculate=False, **kwargs):
        for correspondence in self.utils.reference(pdb):
            ref_pdb = correspondence['pdb1']
            mapping = self.utils.mapping(ref_pdb, pdb)
            ref_loops = self.utils.loops(ref_pdb)
            pdb_loops = self.utils.loops(pdb)
            mapped = self.map_loops(pdb_loops, ref_pdb, mapping)

            if not mapped:
                logger.error("No loops could be mapped")
                continue

            yield self.compare(ref_loops, mapped, correspondence['id'])

    def __call__(self, pdbs, **kwargs):
        if not pdbs:
            raise Exception("No pdbs given")

        for pdb in it.imap(lambda p: p.upper(), pdbs):
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
