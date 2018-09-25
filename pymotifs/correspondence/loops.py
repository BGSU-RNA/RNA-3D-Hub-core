import sys
import traceback

from fr3d import geometry as geo

from pymotifs import core
from pymotifs import utils
from pymotifs import models as mod
from pymotifs.utils.structures import Structure as StructureUtil


class MissingNucleotideException(Exception):
    """Raised when there is a missing nucleotide that we need to map.
    """
    pass


class Loader(core.Loader):
    name = 'correspondence_loops'
    update_gap = False

    def __init__(self, config, maker):
        self._overlaps = {}
        self.utils = StructureUtil(maker)
        self.cif = utils.CifData()
        super(Loader, self).__init__(config, maker)

    def overlap(self, coverage):
        if not self._overlaps:
            with self.session() as session:
                for overlap in session.query(mod.LoopOverlapInfo):
                    self._overlaps[overlap.name] = overlap.id

        if coverage not in self._overlaps:
            raise core.InvalidState("Unknown overlap name: " + coverage)

        return self._overlaps[coverage]

    def nts(self, loop):
        structure = self.cif(loop['pdb']).structure()
        return structure.residues(unit_id=loop['nts'])

    def discrepancy(self, loop1, loop2):
        return geo.discrepancy(self.loops(loop1), self.loops(loop2))

    def map(self, nts, mapping):
        if not mapping:
            raise core.InvalidState("Given empty mapping")

        for nt in nts:
            if nt not in mapping:
                raise core.InvalidState("Missing nt %s" % nt)
            yield mapping[nt]

    def map_loops(self, loops, ref, mapping):
        mapped = []
        for loop in loops:
            try:
                nts = list(self.map(loop['nts'], mapping))
            except MissingNucleotideException as err:
                self.logger.warn("Loop %s cannot be mapped to %s", loop['id'],
                                 ref)
                self.logger.warn("Missing nt %s", str(err))
                continue

            mapped.append({'id': loop['id'], 'nts': nts})

        return mapped

    def coverage(self, loop1, loop2):
        if not loop1:
            raise core.InvalidState("Empty first loop")

        if not loop2:
            raise core.InvalidState("Empty second loop")

        loop_nt1 = set(loop1.get('nts', []))
        if not loop_nt1:
            raise core.InvalidState("First loop had no nts")

        nts = loop2.get('nts', [])
        if not nts:
            raise core.InvalidState("Second loop has no nts")

        mapped = set(loop2['nts'])
        if len(mapped) != len(loop2['nts']):
            raise core.InvalidState("Could not map all loop2 nts")

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
        raise core.InvalidState("This should never occur")

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
                one().\
                id

    def compare_loop(self, ref_loop, loops, corr_id):
        overlapping = []
        for loop in loops:
            try:
                cover = self.coverage(ref_loop, loop)
            except:
                self.logger.error("Exception in overlap between %s, %s",
                                  ref_loop['id'], loop['id'])
                self.logger.error(traceback.format_exc(sys.exc_info()))
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

    def has_data(self, reference, pdb):
        with self.session() as session:
            query = self.__base_info_query__(session, reference, pdb)
            return bool(query.count())

    def remove(self, pdb):
        with self.session() as session:
            session.query(mod.CorrespondenceLoops).\
                join(mod.CorrespondenceInfo,
                     mod.CorrespondenceInfo.correspondence_id == mod.CorrespondenceLoops.correspondence_id).\
                filter(mod.CorrespondenceInfo.pdb_id_1 == pdb).\
                delete()

    def transform(self, pdb):
        cif = self.cif(pdb)
        with self.session() as session:
            query = session.query(mod.CorrespondenceInfo).\
                join(mod.ExpSeqInfo, mod.ExpSeqInfo.exp_seq_id == mod.CorrespondenceInfo.exp_seq_id_1).\
                filter(mod.ExpSeqInfo.pdb_id == pdb)
            return [(cif, result.id) for result in query]

    def data(self, entry, **kwargs):
        # corr_id, loop1, loop2 = entry
        # ref_pdb = correspondence['pdb1']
        # mapping = self.utils.mapping(ref_pdb, pdb)
        # ref_loops = self.utils.loops(ref_pdb)
        # pdb_loops = self.utils.loops(pdb)
        # mapped = self.map_loops(pdb_loops, ref_pdb, mapping)
        # return self.compare(ref_loops, mapped, correspondence['id'])

        return []
