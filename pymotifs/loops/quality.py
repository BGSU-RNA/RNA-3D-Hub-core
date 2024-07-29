"""
Program for importing loop quality assurance data into the RNA 3D Hub database.
Data goes in the table loop_qa
Some loops should be disqualified because they have unresolved or missing
nucleotides. The disqualification codes are:

1 - valid loop
2 - missing nucleotides
3 - modified nucleotides
4 - abnormal chain number
5 - incomplete nucleotides
6 - self-complementary internal loop
7 - Too many symmetry operators
8 - A fictional loop ... one with too high values of RSRZ
9 - A loop with pair outside of the RSRZ range
"""

from collections import defaultdict
from collections import namedtuple as nt

import operator as op
import functools as ft
import collections as coll

from sqlalchemy import asc
from sqlalchemy import desc
from sqlalchemy.orm import aliased

from fr3d.unit_ids import decode

from pymotifs import core
from pymotifs.utils import row2dict
from pymotifs.utils import grouper
from pymotifs import models as mod
from pymotifs.units.incomplete import Entry

from pymotifs.constants import RSRZ_PAIRED_OUTLIERS
from pymotifs.constants import RSRZ_FICTIONAL_CUTOFF

#from pymotifs.loops.release import Loader as ReleaseLoader
from pymotifs.loops.extractor import Loader as InfoLoader
from pymotifs.loops.positions import Loader as PositionLoader
from pymotifs.units.incomplete import Loader as IncompleteLoader
from pymotifs.exp_seq.mapping import Loader as ExpSeqMappingLoader
from pymotifs.exp_seq.positions import Loader as ExpSeqPositionLoader

COMPLEMENTARY = {
    'G': 'C',
    'C': 'G',
    'A': 'U',
    'U': 'A'
}


class AssessmentData(nt('AssessmentData', ['incomplete', 'pairs', 'rsrz'])):
    """Just a value class to store some assessment related data.

    Attributes
    ----------
    incomplete : dict
        Data on what nucleotides are incomplete
    pairs : dict
        Dict mapping from loop id to nts in the loop that have internal
        interactions
    rsrz : dict
        Data on rsrz data for all nucleotides
    """
    pass


class Loader(core.SimpleLoader):
    dependencies = set([InfoLoader, PositionLoader,
                        ExpSeqPositionLoader, ExpSeqMappingLoader,
                        IncompleteLoader])

    allow_no_data = True

    @property
    def table(self):
        return mod.LoopQa

    def to_process(self, pdbs, **kwargs):
        """Convert the list of pdbs to only those PDB's with loops that have not been checked by quality yet. By
        doing this this stage is able assert that data is always produced.

        Parameters
        ----------
        pdbs : list
            List of PDB ids

        Returns
        -------
        pdbs : list
            A list of PDBs from the original list that contain loops and have not been checked for quality yet.
        """

        # accumulate PDB ids for which loops already appear in the loop info table
        with self.session() as session:
            query = session.query(mod.LoopInfo.pdb_id).\
                join(mod.LoopPositions,
                     mod.LoopPositions.loop_id == mod.LoopInfo.loop_id).\
                distinct()
            have_loops = {r.pdb_id for r in query}

        # get list of pdbs with related entries in loop_qa
        with self.session() as session:
            query = session.query(mod.LoopInfo.pdb_id).\
                join(mod.LoopQa,
                     mod.LoopQa.loop_id == mod.LoopInfo.loop_id).\
                distinct()
            checked = {r.pdb_id for r in query}

        #Get list of pdbs with loops
        to_use = sorted(set(pdbs).intersection(have_loops))
        #Get list of pdbs with loops that have NOT been checked for quality yet
        to_use = sorted(set(to_use).difference(checked))

        if not to_use:
            raise core.Skip("All loops in the current PDB list have gone through QA")
        return to_use

    def query(self, session, pdb):
        """Create a query to find all loop qa entries for the given pdb.

        :param session: The session to query with.
        :param str pdb: The PDB id to look up.
        :returns: The query to use.
        """

        return session.query(mod.LoopQa).\
            join(mod.LoopInfo, mod.LoopInfo.loop_id == mod.LoopQa.loop_id).\
            filter(mod.LoopInfo.pdb_id == pdb)

    def paired(self, pdb):
        """Load all within loop basepairs.

        Parameters
        ----------
        pdb : str
            The pdb id to use.

        Returns
        -------
            A dict mapping from loop id to a set of paired units.
        """
        with self.session() as session:
            pos1 = aliased(mod.LoopPositions)
            pos2 = aliased(mod.LoopPositions)
            interactions = mod.UnitPairsInteractions
            info = mod.BpFamilyInfo
            query = session.query(interactions.unit_id_1.label('unit1'),
                                  interactions.unit_id_2.label('unit2'),
                                  pos1.loop_id,
                                  ).\
                join(pos1, pos1.unit_id == interactions.unit_id_1).\
                join(pos2, pos2.unit_id == interactions.unit_id_2).\
                join(info, info.bp_family_id == interactions.f_lwbp).\
                filter(interactions.pdb_id == pdb).\
                filter(interactions.program == 'matlab').\
                filter(pos1.loop_id == pos2.loop_id).\
                filter(info.is_near == 0)

            pairs = defaultdict(set)
            for result in query:
                pairs[result.loop_id].add(result.unit1)
                pairs[result.loop_id].add(result.unit2)
            return pairs

    def incomplete(self, pdb):
        """Load all incomplete nucleotides from the database. This will query
        the unit_incomplete for all incomplete data.

        Parameters
        ----------
        pdb : str
            The pdb id to use.

        Returns
        -------
        incomplete : set
            A set of unit ids that are incomplete.
        """
        with self.session() as session:
            query = session.query(mod.UnitIncomplete.pdb_id,
                                  mod.UnitIncomplete.model,
                                  mod.UnitIncomplete.chain,
                                  mod.UnitIncomplete.number,
                                  mod.UnitIncomplete.unit,
                                  mod.UnitIncomplete.alt_id,
                                  mod.UnitIncomplete.ins_code,
                                  ).\
                filter_by(pdb_id=pdb).\
                filter(mod.UnitIncomplete.model.isnot(None))

        return {Entry(**row2dict(r)) for r in query}

    def loops(self, pdb):
        """Load all the loops in the pdb file.

        Parameters
        ----------
        pdb : str
            The pdb id to use

        Returns
        -------
        loops : list
            A list of loop dictionaries. Each loop will contain an 'nts' list,
            'units' list, chains set, signature string and endpoints list.
        """

        def empty_loop():
            return {
                'nts': [],
                'units': [],
                'chains': set(),
                'signature': [],
                'endpoints': [],
                }

        with self.session() as session:
            badloops = [ 'HL_2H0S_001', 'IL_2H0S_001' ]
            query = session.query(mod.LoopInfo.loop_id.label('id'),
                                  mod.LoopInfo.type,
                                  mod.LoopInfo.pdb_id.label('pdb'),
                                  mod.UnitInfo.unit_id.label('uid'),
                                  mod.UnitInfo.unit.label('unit'),
                                  mod.LoopInfo.seq.label('seq'),
                                  mod.UnitInfo.chain.label('chain'),
                                  mod.UnitInfo.number.label('number'),
                                  mod.LoopPositions.border,
                                  ).\
                join(mod.LoopPositions,
                     mod.LoopPositions.loop_id == mod.LoopInfo.loop_id).\
                join(mod.UnitInfo,
                     mod.UnitInfo.unit_id == mod.LoopPositions.unit_id).\
                filter(mod.LoopInfo.pdb_id == pdb).\
	        	filter(~mod.LoopInfo.loop_id.in_(badloops)).\
                order_by(asc(mod.LoopPositions.position))

            loops = coll.defaultdict(empty_loop)
            for result in query:
                current = loops[result.id]
                current['id'] = result.id
                current['type'] = result.type
                current['nts'].append(result.uid)
                current['units'].append(result.unit)
                current['chains'].add(result.chain)
                current['signature'].append(result.number)
                current['pdb'] = result.pdb
                current['seq'] = result.seq

                if result.border:
                    current['endpoints'].append(result.uid)

                loops[result.id] = current

            loops = loops.values()
            for loop in loops:
                ends = loop['endpoints']
                loop['endpoints'] = [(e1, e2) for (e1, e2) in grouper(2, ends)]
            return sorted(loops, key=op.itemgetter('id'))

    def position_info(self, unit):
        """
        Get the information about a position in an experimental sequence
        using a unit id.
        """

        # self.logger.info("Finding position for %s" % unit)
        try:
            with self.session() as session:
                pos = mod.ExpSeqPosition
                mapping = mod.ExpSeqUnitMapping
                result = session.query(pos.index,
                                       pos.exp_seq_id,
                                       mod.UnitInfo.chain,
                                       mod.UnitInfo.model,
                                       mod.UnitInfo.sym_op,
                                       ).\
                    join(mapping,
                         mapping.exp_seq_position_id == pos.exp_seq_position_id).\
                    join(mod.UnitInfo,
                         mod.UnitInfo.unit_id == mapping.unit_id).\
                    filter(mapping.unit_id == unit).\
                    one()

            return row2dict(result)
        except:
            # handle the case where the unit id in the database table ends with ||A or ||B
            # but that is not being stored in unit.  Not sure why not.
            self.logger.info('Looking up sequence position of alternates of '+unit)
            newunit = '%' + unit + '%'
            with self.session() as session:
                pos = mod.ExpSeqPosition
                mapping = mod.ExpSeqUnitMapping
                result = session.query(pos.index,
                                       pos.exp_seq_id,
                                       mod.UnitInfo.chain,
                                       mod.UnitInfo.model,
                                       mod.UnitInfo.sym_op,
                                       ).\
                    join(mapping,
                         mapping.exp_seq_position_id == pos.exp_seq_position_id).\
                    join(mod.UnitInfo,
                         mod.UnitInfo.unit_id == mapping.unit_id).\
                    filter(mapping.unit_id.like(newunit)).\
                    first()

            if not result:
                self.logger.info('No experimental sequence position for ' + unit)

            return row2dict(result)


    def units_between(self, unit1, unit2):
        """
        Get a list of all units between two units. This assumes they are on
        the same chain and have the same symmetry operator.
        This method refers to the ExpSeqUnitMapping table, which is surprising
        because the same information could be obtained from chain_index
        in the unit_info table.
        """

        start = self.position_info(unit1)
        stop = self.position_info(unit2)

        # self.logger.info('Units are %s and %s' % (unit1, unit2))


        with self.session() as session:
            units = mod.UnitInfo
            mapping = mod.ExpSeqUnitMapping
            pos = mod.ExpSeqPosition
            query = session.query(units.pdb_id,
                                  units.model,
                                  units.chain,
                                  units.number,
                                  units.unit,
                                  units.alt_id,
                                  units.ins_code,
                                  ).\
                join(mapping,
                     mapping.unit_id == units.unit_id).\
                join(pos,
                     mapping.exp_seq_position_id == pos.exp_seq_position_id).\
                filter(pos.exp_seq_id == start['exp_seq_id']).\
                filter(pos.index >= start['index']).\
                filter(pos.index <= stop['index']).\
                filter(units.chain == start['chain']).\
                filter(units.model == start['model']).\
                filter(units.sym_op == start['sym_op']).\
                distinct().\
                order_by(asc(pos.index))

            return [Entry(**row2dict(r)) for r in query]

    def complementary_sequence(self, loop):
        """Detect if a sequence is complementary.
        """

        if len(loop['units']) % 2 != 0:
            return False

        parts = []
        for (start, stop) in loop['endpoints']:
            units = self.units_between(start, stop)
            parts.append([u.unit for u in units])

        pair = zip(parts, reversed(parts))
        for (part1, part2) in pair:
            if len(part1) != len(part2):
                return False
            seq_pair = zip(part1, reversed(part2))
            for (first, second) in seq_pair:
                if second not in COMPLEMENTARY or first != COMPLEMENTARY[second]:
                    return False

        return True

    def has_no_non_cWW(self, loop):
        """
        Check if there are non-cWW interactions within the loop.
        """

        with self.session() as session:
            inters = mod.UnitPairsInteractions
            bps = mod.BpFamilyInfo
            query = session.query(inters).\
                join(bps, bps.bp_family_id == inters.f_lwbp).\
                filter(inters.unit_id_1.in_(loop['nts'])).\
                filter(inters.unit_id_2.in_(loop['nts'])).\
                filter(inters.program == 'python').\
                filter(bps.is_near == 0).\
                filter(bps.bp_family_id != 'cWW')

            return not bool(query.count())

    def is_complementary(self, loop):
        """
        Check if a loop has a complementary sequence. This requires that the
        loop have no non-cWW basepairs (though near are allowed) between them.
        If so then we consider it as a complemenatry loop since it is likely
        just a poorly modeled helix.
        """

        if loop['type'] != 'IL':
            return False

        return self.has_no_non_cWW(loop) and self.complementary_sequence(loop)

    def modified_bases(self, loop):
        """Get a list of all modified bases, if any in this loop.
        """

        modified = []
        normal = ft.partial(op.contains, set(['A', 'C', 'G', 'U']))
        for (start, end) in loop['endpoints']:
            units = self.units_between(start, end)
            modified.extend(u.unit for u in units if not normal(u.unit))

        return modified

    def has_modified(self, loop):
        """Check if there are any modified nucleotides in the loops. These
        cannot yet be processed.
        """
        return bool(self.modified_bases(loop))

    def has_breaks(self, loop):
        """Check if there are any chain breaks within the loop.

        :returns: Bool, true if the loop has any breaks.
        """

        for (u1, u2) in loop['endpoints']:
            start = self.position_info(u1)
            stop = self.position_info(u2)

            with self.session() as session:
                mapping = mod.ExpSeqUnitMapping
                pos = mod.ExpSeqPosition
                query = session.query(mapping).\
                    join(pos,
                         mapping.exp_seq_position_id == pos.exp_seq_position_id).\
                    filter(pos.exp_seq_id == start['exp_seq_id']).\
                    filter(pos.index >= start['index']).\
                    filter(pos.index <= stop['index']).\
                    filter(mapping.chain == start['chain']).\
                    filter(mapping.unit_id == None)

                return bool(query.count())

    def has_incomplete_nucleotides(self, incomplete, loop):
        """Check if any of the nucleotides in the loop are incomplete, that is
        are missing atoms that they should not be.
        """

        for (u1, u2) in loop['endpoints']:
            for unit in self.units_between(u1, u2):
                if unit in incomplete:
                    return True
        return False

    def bad_chain_number(self, loop):
        """Check if there are the correct number of chains for each loop.
        """

        if not loop['chains']:
            return True

        if loop['type'] == 'HL':
            return len(loop['chains']) != 1
        if loop['type'] == 'IL':
            return len(loop['chains']) > 2
        if loop['type'] == 'J3':
            return len(loop['chains']) > 3
        raise core.InvalidState("Unknown type of loop")

    def too_many_sym_ops(self, loop):
        """Detect if we have only 1 symmetry operator in this loop. This is
        an invalid loop.

        :param dict loop: The loop to examine.
        :returns: Bool, true if there is more symmetry operator.
        """

        ops = set()
        for unit in loop['nts']:
            parts = decode(unit)
            ops.add(parts['symmetry'])
        if loop['type'] == 'HL':
            return len(ops) != 1
        if loop['type'] == 'IL':
            return len(ops) > 2
        if loop['type'] == 'J3':
            return len(ops) > 3

    def rsrz_data(self, pdb):
        """Load the unit quality data for all units in the given structure.

        Parameters
        ----------
        pdb : str
            The PDB id to use.

        Returns
        -------
        data : dict
            A dict mapping from unit id to the rsr z_score.
            Removed the following line because it crashed:
#                filter_by(unit.in_('A','C','G','U')).\
        """
        data = defaultdict(lambda: None)
        na_types = ['rna','dna','hybrid']
        with self.session() as session:
            query = session.query(mod.UnitQuality).\
                join(mod.UnitInfo,
                     mod.UnitInfo.unit_id == mod.UnitQuality.unit_id).\
                filter(mod.UnitInfo.unit_type_id.in_(na_types)).\
                filter_by(pdb_id=pdb)
            data.update({r.unit_id: r.real_space_r_z_score for r in query})
            return data

    def is_fictional_loop(self, rsrz, loop):
        """
        Detect if this loop is fictional. A fictional loop is defined as one
        where all nucleotides in the loop have an RSRZ value are above some cutoff.
        The cutoff used here is the RSRZ_FICTIONAL_CUTOFF imported from
        pymotifs.constants.

        Parameters
        ----------
        rsrz : dict
            A dict mapping from unit id to rsrz score.
        loop : dict
            A loop dict with a nt entries that contains the list of all unit
            ids in the loop.

        Returns
        -------
        is_fictional : bool
            True if all nucleotides in the loop have an RSRZ greater than the
            cutoff.
        """

        return all(rsrz[unit] is None or rsrz[unit] >= RSRZ_FICTIONAL_CUTOFF for unit in loop['nts'])

    def is_fictional_pair(self, pairs, rsrz, loop):
        """
        """
        paired = pairs[loop['id']]
        return any(rsrz[u] >= RSRZ_PAIRED_OUTLIERS for u in loop['nts'] if u in paired)

    def status(self, assess, loop):
        """Compute the status code. The status code is defined above.

        Parameters
        ----------
        assess : AssessmentData
            An assessment data object for information about the structure
        loop : dict
            A dictionary for containing information about the loop.

        Returns
        -------
        status : int
            The status code
        """

        if self.has_breaks(loop):
            return 2
        if self.has_modified(loop):
            return 3
        if self.bad_chain_number(loop):
            return 4
        if self.too_many_sym_ops(loop):
            return 7
        if self.has_incomplete_nucleotides(assess.incomplete, loop):
            return 5
        if self.is_complementary(loop):
            return 6
        if self.is_fictional_loop(assess.rsrz, loop):
            return 8
# The following check is not able to work; need to pass in pairs, how to get it?
#        if self.is_fictional_pair(assess.rsrz, loop):
#            return 9
        return 1

    def quality(self, assess, loop):
        """Compute the quality information for the given loop.

        Parameters
        ----------
        assess : AssessmentData
            The assessment data to use.
        loop : dict
            A loop dictionary to use.

        Returns
        -------
        quality : dict
            A quality dict to write to the database.
        """

        self.logger.debug("Examining loop %s", str(loop))
        seq = None
        if self.is_complementary(loop):
            seq = loop['seq'].replace('*', ',')

        mods = None
        if self.modified_bases(loop):
            mods = ', '.join(self.modified_bases(loop))

        return {
            'loop_id': loop['id'],
            'status': self.status(assess, loop),
            'modifications': mods,
            'nt_signature': ', '.join(str(s) for s in loop['signature']),
            'complementary': seq
        }

    def assessment_data(self, pdb):
        return AssessmentData(incomplete=self.incomplete(pdb),
                              pairs=self.paired(pdb),
                              rsrz=self.rsrz_data(pdb))

    def data(self, pdb, **kwargs):
        """Compute the qa status of each loop in the structure.

        :param str pdb: The pdb id to use.
        :returns: A list of the status entries for all loops in the structure.
        """

        assess = self.assessment_data(pdb)
        return [self.quality(assess, l) for l in self.loops(pdb)]
