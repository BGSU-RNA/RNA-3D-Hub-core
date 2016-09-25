"""
Program for importing loop quality assurance data into the RNA 3D Hub database.

Some loops should be disqualified because they have unresolved or missing
nucleotides. The disqualification codes are:

1 - valid loop
2 - missing nucleotides
3 - modified nucleotides
4 - abnormal chain number
5 - incomplete nucleotides
6 - self-complementary internal loop
7 - Too many symmetry operators
"""

import operator as op
import functools as ft
import itertools as it
import collections as coll

from sqlalchemy import asc
from sqlalchemy import desc

from fr3d.unit_ids import decode

from pymotifs import core
from pymotifs.utils import row2dict
from pymotifs.utils import grouper
from pymotifs import models as mod

from pymotifs.loops.release import Loader as ReleaseLoader
from pymotifs.loops.extractor import Loader as InfoLoader
from pymotifs.loops.positions import Loader as PositionLoader
from pymotifs.exp_seq.mapping import Loader as ExpSeqMappingLoader
from pymotifs.exp_seq.positions import Loader as ExpSeqPositionLoader

COMPLEMENTARY = {
    'G': 'C',
    'C': 'G',
    'A': 'U',
    'U': 'A'
}


class Loader(core.SimpleLoader):
    dependencies = set([ReleaseLoader, InfoLoader, PositionLoader,
                        ExpSeqPositionLoader, ExpSeqMappingLoader])

    @property
    def table(self):
        return mod.LoopQa

    def to_process(self, pdbs, **kwargs):
        """Convert the list of pdbs to only those PDB's that have a loop. By
        doing this this stage is able assert that this always data produced.

        Parameters
        ----------
        pdbs : list
            List of PDB ids

        Returns
        -------
        pdbs : list
            A list of PDBs from the original list that contain loops.
        """

        with self.session() as session:
            query = session.query(mod.LoopInfo.pdb_id).\
                distinct()
            known = set(r.pdb_id for r in query)

        return sorted(set(pdbs).intersection(known))

    def current_id(self):
        """Compute the current loop release id.

        :returns: The current loop release id string, eg 1.2.
        """
        with self.session() as session:
            query = session.query(mod.LoopReleases.loop_release_id).\
                order_by(desc(mod.LoopReleases.date)).\
                limit(1)

            if not query.count():
                return '0.0'

            return query.one().loop_release_id

    def query(self, session, pdb):
        """Create a query to find all loop qa entries for the given pdb.

        :param session: The session to query with.
        :param str pdb: The PDB id to look up.
        :returns: The query to use.
        """

        release_id = self.current_id()
        return session.query(mod.LoopQa).\
            join(mod.LoopInfo, mod.LoopInfo.loop_id == mod.LoopQa.loop_id).\
            filter(mod.LoopInfo.pdb_id == pdb).\
            filter(mod.LoopQa.loop_release_id == release_id)

    def components(self, pdb):
        structure = self.structure(pdb)
        components = {}
        for comp in structure.residues():
            components[comp.unit_id()] = comp
        return components

    def incomplete(self, cif):
        data = set()
        total = it.chain(getattr(cif, 'pdbx_unobs_or_zero_occ_residues', []),
                         getattr(cif, 'pdbx_unobs_or_zero_occ_atoms', []))
        for row in total:
            model = int(row['PDB_model_num'])
            chain = row['auth_asym_id']
            num = int(row['auth_seq_id'])
            seq = row['auth_comp_id']
            ins = None
            if row['PDB_ins_code'] != '?':
                ins = row['PDB_ins_code']

            alt_id = None
            if 'label_alt_id' in row and row['label_alt_id'] != '?':
                alt_id = row['label_alt_id']

            key = (model, chain, num, seq, alt_id, ins)
            data.add(key)

        return data

    def loops(self, pdb):
        """Load all the loops in the pdb file.
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
            return loops

    def position_info(self, unit):
        """Get the information about a position in an experimental sequence
        using a unit id.
        """

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

    def units_between(self, unit1, unit2):
        """Get a list of all units between two units. This assumes they are on
        the same chain and have the same symmetry operator.
        """

        start = self.position_info(unit1)
        stop = self.position_info(unit2)
        with self.session() as session:
            units = mod.UnitInfo
            mapping = mod.ExpSeqUnitMapping
            pos = mod.ExpSeqPosition
            query = session.query(units).\
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

            return [row2dict(r) for r in query]

    def complementary_sequence(self, loop):
        """Detect if a sequence is complementary.
        """

        if len(loop['units']) % 2 != 0:
            return False

        parts = []
        for (start, stop) in loop['endpoints']:
            units = self.units_between(start, stop)
            parts.append([u['unit'] for u in units])

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
        """Check if there are non-cWW interactions within the loop.
        """

        with self.session() as session:
            inters = mod.UnitPairsInteractions
            bps = mod.BpFamilyInfo
            query = session.query(inters).\
                join(bps, bps.bp_family_id == inters.f_lwbp).\
                filter(inters.unit_id_1.in_(loop['nts'])).\
                filter(inters.unit_id_2.in_(loop['nts'])).\
                filter(bps.is_near == 0).\
                filter(bps.bp_family_id != 'cWW')

            return not bool(query.count())

    def is_complementary(self, loop):
        """Check if a loop has a complementary sequence. This requires that the
        loop have no non-cWW basepairs (though near are allowed) between them.
        If so then we consider it as a complemenatry loop since it is likely
        just a poorly modeled helix.
        """

        return self.has_no_non_cWW(loop) and self.complementary_sequence(loop)

    def modified_bases(self, loop):
        """Get a list of all modified bases, if any in this loop.
        """

        modified = []
        normal = ft.partial(op.contains, set(['A', 'C', 'G', 'U']))
        for (start, end) in loop['endpoints']:
            units = self.units_between(start, end)
            modified.extend(u['unit'] for u in units if not normal(u['unit']))

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
                key = (int(unit['model']), unit['chain'], int(unit['number']),
                       unit['unit'], unit['ins_code'], unit['alt_id'])
                if key in incomplete:
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
        for nt in loop['nts']:
            parts = decode(nt)
            ops.add(parts['symmetry'])
        return len(ops) != 1

    def status(self, incomplete, loop):
        """Compute the status code. The status code is defined above.

        :param dict loop: A loop as returned from loops.
        :returns: The status code.
        """

        if self.has_breaks(loop):
            return 2
        if self.has_modified(loop):
            return 3
        if self.bad_chain_number(loop):
            return 4
        if self.too_many_sym_ops(loop):
            return 7
        if self.has_incomplete_nucleotides(incomplete, loop):
            return 5
        if self.is_complementary(loop):
            return 6
        return 1

    def quality(self, incomplete, release_id, loop):
        seq = None
        if self.is_complementary(loop):
            seq = loop['seq'].replace('*', ',')

        mods = None
        if self.modified_bases(loop):
            mods = ', '.join(self.modified_bases(loop))

        return {
            'loop_id': loop['id'],
            'status': self.status(incomplete, loop),
            'modifications': mods,
            'nt_signature': ', '.join(str(s) for s in loop['signature']),
            'complementary': seq,
            'loop_release_id': release_id,
        }

    def data(self, pdb, **kwargs):
        """Compute the qa status of each loop in the structure.

        :param str pdb: The pdb id to use.
        :returns: A list of the status entries for all loops in the structure.
        """

        cif = self.cif(pdb)
        partial = self.incomplete(cif)
        release_id = self.current_id()
        return [self.quality(partial, release_id, l) for l in self.loops(pdb)]
