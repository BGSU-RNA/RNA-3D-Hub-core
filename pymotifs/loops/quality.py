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
from fr3d.modified.mapping import modified_base_to_parent

from pymotifs import core
from pymotifs.utils import row2dict
from pymotifs.utils import grouper
from pymotifs import models as mod
from pymotifs.units.incomplete import Entry

from pymotifs.constants import RSRZ_PAIRED_OUTLIERS
from pymotifs.constants import RSRZ_FICTIONAL_CUTOFF

#from pymotifs.loops.release import Loader as ReleaseLoader
from pymotifs.loops.save_loops import Loader as SaveLoopsLoader
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
    """
    Just a value class to store some assessment related data.

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
    dependencies = set([SaveLoopsLoader,
                        ExpSeqPositionLoader, ExpSeqMappingLoader,
                        IncompleteLoader])

    allow_no_data = True
    mark = False
    use_marks = False

    @property
    def table(self):
        return mod.LoopQa

    def to_process(self, pdbs, **kwargs):
        """
        Convert the list of pdbs to only those PDB's with loops that have not been checked by quality yet.

        Parameters
        ----------
        pdbs : list
            List of PDB ids

        Returns
        -------
        pdbs : list
            A list of PDBs from the original list that contain loops and have not been checked for quality yet.
        """

        # accumulate loop ids in the loop_info and loop_positions from pdb_id in pdbs
        with self.session() as session:
            query = session.query(mod.LoopInfo.loop_id).\
                                join(mod.LoopPositions,
                                    mod.LoopPositions.loop_id == mod.LoopInfo.loop_id).\
                                filter(mod.LoopInfo.pdb_id.in_(pdbs)).\
                                distinct()
            all_loops = {r.loop_id for r in query}

        print('all_loops:', all_loops)

        # get list of loop_ids in loop_qa from pdbs
        with self.session() as session:
            query = session.query(mod.LoopInfo.loop_id).\
                            join(mod.LoopQa,
                                mod.LoopQa.loop_id == mod.LoopInfo.loop_id).\
                            filter(mod.LoopInfo.pdb_id.in_(pdbs)).\
                            distinct()
            checked_loops = {r.loop_id for r in query}

        print('checked_loops already in loop_qa:', checked_loops)

        # get list of loops that have not been checked
        not_checked_loops = all_loops - checked_loops

        # get list of pdbs that have at least one loop that has not been checked
        to_use = set()
        for loop_id in not_checked_loops:
            fields = loop_id.split('_')
            to_use.add(fields[1])

        if not to_use:
            raise core.Skip("All loops in the current PDB list have gone through QA")
        return sorted(to_use)

    def query(self, session, pdb):
        """
        Create a query to find all loop qa entries for the given pdb.

        :param session: The session to query with.
        :param str pdb: The PDB id to look up.
        :returns: The query to use.
        """

        # since to_process correctly identifies the pdb ids to use,
        # make sure to pass back an empty query so it thinks work needs to be done
        return session.query(mod.LoopQa).filter(mod.LoopQa.loop_id == 'aaa')

        return session.query(mod.LoopQa).\
                            join(mod.LoopInfo, mod.LoopInfo.loop_id == mod.LoopQa.loop_id).\
                            filter(mod.LoopInfo.pdb_id == pdb)

    def paired(self, pdb):
        """
        Load all within loop basepairs.

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
            interactions = mod.UnitPairsInteractions2024
            info = mod.BpFamilyInfo
            query = session.query(interactions.unit_id_1.label('unit1'),
                                  interactions.unit_id_2.label('unit2'),
                                  pos1.loop_id,
                                  ).\
                join(pos1, pos1.unit_id == interactions.unit_id_1).\
                join(pos2, pos2.unit_id == interactions.unit_id_2).\
                join(info, info.bp_family_id == interactions.f_lwbp).\
                filter(interactions.pdb_id == pdb).\
                filter(interactions.program == 'fr3d').\
                filter(pos1.loop_id == pos2.loop_id).\
                filter(info.is_near == 0)

            pairs = defaultdict(set)
            for result in query:
                pairs[result.loop_id].add(result.unit1)
                pairs[result.loop_id].add(result.unit2)
            return pairs

    def incomplete(self, pdb):
        """
        Load all incomplete nucleotides from the database. This will query
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
        """
        Load all the loops in the pdb file that are not already in loop_qa table.

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
                'position_2023': [],
                'chain_index': [],
                'border': [],
                }

        # get list of loop ids from this pdb already in loop_qa
        with self.session() as session:
            query = session.query(mod.LoopInfo.loop_id).\
                join(mod.LoopQa,
                     mod.LoopQa.loop_id == mod.LoopInfo.loop_id).\
                filter(mod.LoopInfo.pdb_id == pdb).\
                distinct()
            checked = {r.loop_id for r in query}

        with self.session() as session:
            badloops = [ 'HL_2H0S_001', 'IL_2H0S_001' ]
            query = session.query(mod.LoopInfo.loop_id.label('id'),
                                  mod.LoopInfo.type,
                                  mod.LoopInfo.pdb_id.label('pdb'),
                                  mod.UnitInfo.unit_id.label('uid'),
                                  mod.UnitInfo.unit.label('unit'),
                                  mod.UnitInfo.chain.label('chain'),
                                  mod.UnitInfo.number.label('number'),
                                  mod.UnitInfo.chain_index,
                                  mod.LoopInfo.seq.label('seq'),
                                  mod.LoopPositions.border,
                                  mod.LoopPositions.position_2023
                                  ).\
                join(mod.LoopPositions,
                     mod.LoopPositions.loop_id == mod.LoopInfo.loop_id).\
                join(mod.UnitInfo,
                     mod.UnitInfo.unit_id == mod.LoopPositions.unit_id).\
                filter(mod.LoopInfo.pdb_id == pdb).\
	        	filter(~mod.LoopInfo.loop_id.in_(badloops)).\
                filter(~mod.LoopInfo.loop_id.in_(checked)).\
                order_by(asc(mod.LoopPositions.position_2023))

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
                current['position_2023'].append(result.position_2023)
                current['chain_index'].append(result.chain_index)
                current['border'].append(result.border)

                # if result.border:
                #     current['endpoints'].append(result.uid)

                loops[result.id] = current

            loops = loops.values()

            # for loop in loops:
            #     ends = loop['endpoints']
            #     if len(ends) % 2 != 0:
            #         loop['endpoints'] = []
            #     else:
            #         loop['endpoints'] = [(e1, e2) for (e1, e2) in grouper(2, ends)]

            # loops = [l for l in loops if l['endpoints']]

            return sorted(loops, key=op.itemgetter('id'))

    def position_info(self, unit):
        """
        Get the information about a position in an experimental sequence
        using a unit id.
        """

        self.logger.info("Finding position for %s" % unit)
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


    def unit_info(self, unit_id):
        """
        Get the information about a unit using a unit id.
        """

        with self.session() as session:
            units = mod.UnitInfo
            query = session.query(units.pdb_id,
                                  units.model,
                                  units.chain,
                                  units.number,
                                  units.unit,
                                  units.alt_id,
                                  units.ins_code,
                                  units.chain_index,
                                  units.sym_op
                                  ).\
                filter(units.unit_id == unit_id)

            return row2dict(query.one())


    def units_between(self, unit1, unit2):
        """
        Get a list of all units between two units. This assumes they are on
        the same chain and have the same symmetry operator.
        Use chain_index from unit_info.
        """

        # plan: split the fields, get a range of integers, get unit ids satisfying those
        # that's more robust

        unit1_info = self.unit_info(unit1)
        unit2_info = self.unit_info(unit2)

        if unit1_info["pdb_id"] != unit2_info["pdb_id"]:
            return []
        if unit1_info["model"] != unit2_info["model"]:
            return []
        if unit1_info["chain"] != unit2_info["chain"]:
            return []
        if unit1_info["sym_op"] != unit2_info["sym_op"]:
            return []

        start = int(unit1_info["chain_index"])
        stop  = int(unit2_info["chain_index"])

        self.logger.info('Units are %s and %s' % (unit1, unit2))
        self.logger.info('Chain index range is %d to %d' % (start, stop))

        between = list(range(start+1, stop))

        with self.session() as session:
            units = mod.UnitInfo
            query = session.query(units.pdb_id,
                                  units.model,
                                  units.chain,
                                  units.number,
                                  units.unit,
                                  units.alt_id,
                                  units.ins_code,
                                  ).\
                filter(units.pdb_id == unit1_info["pdb_id"]).\
                filter(units.model == unit1_info["model"]).\
                filter(units.chain == unit1_info["chain"]).\
                filter(units.sym_op == unit1_info["sym_op"]).\
                filter(units.chain_index.in_(between)).\
                distinct().\
                order_by(asc(units.chain_index))

            return [Entry(**row2dict(r)) for r in query]


    # def complementary_sequence(self, loop):
    #     """
    #     Detect if an IL sequence is complementary.

    #     Trouble:  methodology fails on CUC,GCGAG because its reversed version
    #     is GAGCG,CUC and the entire sequence is complementary, if you ignore
    #     the commas.

    #     New python code does not extract complementary IL, so skip this.
    #     """

    #     # must have an even number of nucleotides to proceed
    #     if len(loop['units']) % 2 != 0:
    #         return False

    #     parts = []
    #     for (start, stop) in loop['endpoints']:
    #         units = self.units_between(start, stop)
    #         for unit in units:
    #             u = unit.unit
    #             if u in ['A','C','G','U']:
    #                 parts.append(u)
    #             elif u in ['DA','DC','DG']:
    #                 parts.append(u[1])
    #             elif u in ['DT']:
    #                 parts.append('U')
    #             elif u in modified_base_to_parent:
    #                 # loop up parent nucleotide
    #                 p = modified_base_to_parent[u]
    #                 if p in ['A','C','G','U']:
    #                     parts.append(p)
    #                 elif p in ['DA','DC','DG']:
    #                     parts.append(p[1])
    #                 elif p in ['DT']:
    #                     parts.append('U')
    #                 else:
    #                     parts.append(u)
    #             else:
    #                 parts.append(u)

    #     pair = zip(parts, reversed(parts))
    #     for (part1, part2) in pair:
    #         if len(part1) != len(part2):
    #             return False
    #         seq_pair = zip(part1, reversed(part2))
    #         for (first, second) in seq_pair:
    #             if second not in COMPLEMENTARY or first != COMPLEMENTARY[second]:
    #                 return False

    #     return True

    def has_no_non_cWW(self, loop):
        """
        Check if there are non-cWW interactions within the loop.
        """

        with self.session() as session:
            inters = mod.UnitPairsInteractions2024
            bps = mod.BpFamilyInfo
            query = session.query(inters).\
                join(bps, bps.bp_family_id == inters.f_lwbp).\
                filter(inters.unit_id_1.in_(loop['nts'])).\
                filter(inters.unit_id_2.in_(loop['nts'])).\
                filter(inters.program == 'fr3d').\
                filter(bps.is_near == 0).\
                filter(bps.bp_family_id != 'cWW')

            return not bool(query.count())

    # def is_complementary(self, loop):
    #     """
    #     Check if a loop has a complementary sequence. This requires that the
    #     loop have no non-cWW basepairs (though near are allowed) between them.
    #     If so then we consider it as a complemenatry loop since it is likely
    #     just a poorly modeled helix.
    #     """

    #     if loop['type'] != 'IL':
    #         return False

    #     return self.has_no_non_cWW(loop) and self.complementary_sequence(loop)

    def modified_bases(self, loop):
        """
        Get a list of all modified bases, if any, in this loop.
        """

        modified = []

        # normal = ft.partial(op.contains, set(['A', 'C', 'G', 'U','DA','DC','DG','DT']))
        # for (start, end) in loop['endpoints']:
        #     units = self.units_between(start, end)
        #     modified.extend(u.unit for u in units if not normal(u.unit))

        for unit in loop['units']:
            if not unit in ['A','C','G','U','DA','DC','DG','DT']:
                modified.append(unit)

        return modified

    def has_modified(self, loop):
        """
        Check if there are any modified nucleotides in the loops.
        """
        return bool(self.modified_bases(loop))

    # def has_breaks_old(self, loop):
    #     """
    #     Check if there are any chain breaks within the loop.

    #     :returns: Bool, true if the loop has any breaks.
    #     """

    #     for (u1, u2) in loop['endpoints']:
    #         start = self.position_info(u1)
    #         stop = self.position_info(u2)

    #         with self.session() as session:
    #             mapping = mod.ExpSeqUnitMapping
    #             pos = mod.ExpSeqPosition
    #             query = session.query(mapping).\
    #                 join(pos,
    #                      mapping.exp_seq_position_id == pos.exp_seq_position_id).\
    #                 filter(pos.exp_seq_id == start['exp_seq_id']).\
    #                 filter(pos.index >= start['index']).\
    #                 filter(pos.index <= stop['index']).\
    #                 filter(mapping.chain == start['chain']).\
    #                 filter(mapping.unit_id == None)

    #             return bool(query.count())


    def has_breaks_also_old(self, loop):
        """
        Check if there are any chain breaks within the loop.

        :returns: Bool, true if the loop has any breaks.
        """

        for (u1, u2) in loop['endpoints']:
            unit1_info = self.unit_info(u1)
            unit2_info = self.unit_info(u2)

            between = list(range(int(unit1_info["chain_index"])+1, int(unit2_info["chain_index"])))

            with self.session() as session:
                units = mod.UnitInfo
                query = session.query(units.chain_index).\
                    filter(units.pdb_id == unit1_info["pdb_id"]).\
                    filter(units.model == unit1_info["model"]).\
                    filter(units.chain == unit1_info["chain"]).\
                    filter(units.sym_op == unit1_info["sym_op"]).\
                    filter(units.chain_index.in_(between)).\
                    distinct().\
                    order_by(asc(units.chain_index))

                # compute maximum difference between chain_index values
                # if it is greater than 1, then there is a break
                chain_indices = [r.chain_index for r in query]
                for i in range(1,len(chain_indices)):
                    if chain_indices[i] - chain_indices[i-1] > 1:
                        return True

        return False


    def has_breaks(self, loop):
        """
        Check if there are any chain breaks within the loop.

        :returns: Bool, true if the loop has any breaks.
        """

        border = 0
        for i in range(len(loop['nts'])):
            chain_index = loop['chain_index'][i]
            if border % 2 == 1:
                if chain_index - previous_chain_index > 1:
                    return True

            border += loop['border'][i]
            previous_chain_index = chain_index

        return False

    def has_incomplete_nucleotides_old(self, incomplete, loop):
        """
        Check if any of the nucleotides in the loop are incomplete, that is
        are missing atoms that they should not be.
        """

        for (u1, u2) in loop['endpoints']:
            for unit in self.units_between(u1, u2):
                if unit in incomplete:
                    return True
        return False

    def has_incomplete_nucleotides(self, incomplete, loop):
        """
        Check if any of the nucleotides in the loop are incomplete, that is
        are missing atoms that they should not be.
        """

        if len(set(loop['nts']) & incomplete) > 0:
            return True
        return False

    def bad_chain_number(self, loop):
        """
        Check if there are the correct number of distinct chains for each loop.
        For example, if an IL has nucleotides from 3 or more different chains.
        It is not counting strands.
        """

        if not loop['chains']:
            return True

        if loop['type'] == 'HL':
            return len(loop['chains']) > 1
        if loop['type'] == 'IL':
            return len(loop['chains']) > 2
        if loop['type'][0] == 'J':
            return len(loop['chains']) > int(loop['type'].replace("J",""))

        raise core.InvalidState("Unknown type of loop")

    def too_many_sym_ops(self, loop):
        """
        Detect if we have more symmetry operators than chains.
        This is an invalid loop.

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
        if loop['type'][0] == 'J':
            return len(ops) > int(loop['type'].replace("J",""))

    def rsrz_data(self, pdb):
        """
        Load the unit quality data for all units in the given structure.

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

        found_rsrz_value = False

        data = defaultdict(lambda: None)
        na_types = ['rna','dna','hybrid']
        with self.session() as session:
            query = session.query(mod.UnitQuality).\
                join(mod.UnitInfo,
                     mod.UnitInfo.unit_id == mod.UnitQuality.unit_id).\
                filter(mod.UnitInfo.unit_type_id.in_(na_types)).\
                filter_by(pdb_id=pdb)

            for r in query:
                if r.real_space_r_z_score is not None:
                    found_rsrz_value = True
                data[r.unit_id] = r.real_space_r_z_score

            if found_rsrz_value:
                return data
            else:
                return None


    def all_high_RSRZ(self, rsrz, loop):
        """
        Detect if all nucleotides in this loop have RSRZ values that are above some cutoff.
        The cutoff used here is the RSRZ_FICTIONAL_CUTOFF imported from
        pymotifs.constants.
        Loops like that do not have enough experimental basis to consider.

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

        # self.logger.info(rsrz)

        if rsrz:
            return all(rsrz[unit] is None or rsrz[unit] >= RSRZ_FICTIONAL_CUTOFF for unit in loop['nts'])
        else:
            # no rsrz values for this structure, like with cryo-EM and NMR
            return False

    def is_fictional_pair(self, pairs, rsrz, loop):
        """
        """

        if rsrz:
            paired = pairs[loop['id']]
            return any(rsrz[u] >= RSRZ_PAIRED_OUTLIERS for u in loop['nts'] if u in paired)
        else:
            return False

    def status(self, assess, loop):
        """
        Compute the status code. The status code is defined above.

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

        # put the more serious ones first
        if self.has_breaks(loop):
            return 2
        if self.bad_chain_number(loop):
            return 4
        if self.too_many_sym_ops(loop):
            return 7
        if self.has_incomplete_nucleotides(assess.incomplete, loop):
            return 5
        # if self.is_complementary(loop):
        #     # see comment in complementary_sequence about why to skip this
        #     return 6
        if self.all_high_RSRZ(assess.rsrz, loop):
            return 8
        if self.has_modified(loop):
            return 3
        # The following check is not able to work; need to pass in pairs, how to get it?
        # if self.is_fictional_pair(assess.rsrz, loop):
        #     return 9
        return 1

    def quality(self, assess, loop):
        """
        Compute the quality information for the given loop.

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
        # if self.is_complementary(loop):
        #     seq = loop['seq'].replace('*', ',')

        mb = self.modified_bases(loop)
        if mb:
            mods = ', '.join(mb)
        else:
            mods = None

        return {
            'loop_id': loop['id'],
            'status': self.status(assess, loop),
            'modifications': mods,
            # 'nt_signature': ', '.join(str(s) for s in loop['signature']),
            # 'complementary': seq
        }

    def assessment_data(self, pdb):
        return AssessmentData(incomplete=self.incomplete(pdb),
                              pairs=self.paired(pdb),
                              rsrz=self.rsrz_data(pdb))

    def data(self, pdb, **kwargs):
        """
        Compute the qa status of each loop in the structure.

        :param str pdb: The pdb id to use.
        :returns: A list of the status entries for all loops in the structure.
        """

        assess = self.assessment_data(pdb)
        return [self.quality(assess, l) for l in self.loops(pdb)]
