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
7 - cyro-em structure
"""

import collections as coll

from sqlalchemy import asc
from sqlalchemy import desc

from pymotifs import core
from pymotifs import models as mod

from pymotifs.loops.release import Loader as ReleaseLoader
from pymotifs.loops.extractor import Loader as InfoLoader
from pymotifs.loops.positions import Loader as PositionLoader
from pymotifs.exp_seq.mapping import Loader as ExpSeqMappingLoader

COMPLEMENTARY = {
    'G': 'C',
    'C': 'G',
    'A': 'U',
    'U': 'A'
}


class Loader(core.SimpleLoader):
    dependencies = set([ReleaseLoader, InfoLoader, PositionLoader,
                        ExpSeqMappingLoader])

    def current_id(self):
        """Compute the current loop release id.

        :returns: The current loop release id string, eg 1.2.
        """
        with self.session() as session:
            query = session.query(mod.LoopReleases.loop_release_id).\
                order_by(desc(mod.LoopReleases.date)).\
                limit(1)

            if query.count() == 0:
                return '0.0'

            return query.one().loop_releases_id

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

    def loops(self, pdb):
        """Load all the loops in the pdb file.
        """

        def empty_loop():
            return {
                'nts': [],
                'positions': [],
                'chains': set(),
                'signature': [],
                'endpoints': []
            }

        with self.session() as session:
            query = session.query(mod.LoopInfo.loop_id.label('id'),
                                  mod.LoopInfo.type,
                                  mod.LoopInfo.pdb_id.label('pdb'),
                                  mod.UnitInfo.unit_id.label('uid'),
                                  mod.UnitInfo.unit.label('unit'),
                                  mod.LoopInfo.seq.label('seq'),
                                  mod.UnitInfo.chain_name.label('chain'),
                                  mod.UnitInfo.number.label('number'),
                                  mod.LoopPositions.flanking,
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

                if result.flanking:
                    current['endpoints'].append('uid')

                loops[result.id] = current
            return loops.values()

    def complementary_sequence(self, loop):
        """Detect if a sequence is complementary.
        """

        sequence = loop['units']

        midpoint = None
        if len(loop['units']) % 2 != 0:
            midpoint = int(len(loop['units']) / 2)

        pair = zip(sequence, reversed(sequence))
        for index, (first, second) in enumerate(pair):
            if midpoint is not None and index == midpoint:
                continue
            if first != COMPLEMENTARY[second]:
                return False

        return True

    def has_no_non_cWW(self, loop):
        """Check if there are non-cWW interactions within the loop.
        """

        with self.session() as session:
            inters = mod.UnitPairInteractions
            bps = mod.BasePairFamilyInfo
            query = session.query(inters).\
                join(bps, bps.bp_family_id == inters.f_lwbp).\
                filter(inters.unit_id_1.in_(loop['nts'])).\
                filter(inters.unit_id_2.in_(loop['nts'])).\
                filter(bps.is_near == 0).\
                filter(bps.bp_family_id != 'cWW')

            return bool(query.count())

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
        def is_modified(unit):
            return unit not in set(['A', 'C', 'G', 'U'])

        return [unit for unit in loop['units'] if is_modified(unit)]

    def has_modified(self, loop):
        """Check if there are any modified nucleotides in the loops. These
        cannot yet be processed.
        """
        return bool(self.modified_bases(loop))

    def has_breaks(self, loop):
        """Check if there are any chain breaks within the loop.
        """
        pass

    def has_incomplete_nucleotides(self, cif, loop):
        """Check if any of the nucleotides in the loop are incomplete, that is
        are missing atoms that they should not be.
        """
        pass

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

    def is_cyro(self, loop):
        pass

    def status(self, cif, loop):
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
        if self.has_incomplete_nucleotides(cif, loop):
            return 5
        if self.is_complementary(loop):
            return 6
        if self.is_cyro(loop):
            return 7
        return 1

    def data(self, pdb):
        """Compute the qa status of each loop in the structure.

        :param str pdb: The pdb id to use.
        :returns: A list of the status entries for all loops in the structure.
        """

        data = []
        cif = self.cif(pdb)
        release_id = self.current_id()
        for loop in self.loops(pdb):
            data.append({
                'loop_id': loop['id'],
                'status': self.status(cif, loop),
                'modifications': ', '.join(self.modified_bases(loop)),
                'nt_signature': ', '.join(loop['signature']),
                'complementary': self.is_complementary(loop),
                'loop_release_id': release_id,
            })

        return data
