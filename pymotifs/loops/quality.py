"""Program for importing loop quality assurance data into the RNA 3D Hub
database.

Some loops should be disqualified because they have unresolved or missing
nucleotides.
"""

import itertools as it

# from sqlalchemy.orm import aliased

from pymotifs import core
from pymotifs import utils as ut
from pymotifs import models as mod

from pymotifs.loops.extractor import Loader as InfoLoader
from pymotifs.loops.positions import Loader as PositionLoader
from pymotifs.exp_seq.positions import Loader as ExpSeqPositionLoader
from pymotifs.exp_seq.mapping import Loader as ExpSeqUnitLoader
from pymotifs.units.info import Loader as UnitLoader


class Loader(core.SimpleLoader):
    dependencies = set([InfoLoader, ExpSeqPositionLoader, ExpSeqUnitLoader,
                        UnitLoader, PositionLoader])
    table = mod.LoopQa

    complementary = {
        'A': 'U',
        'C': 'G',
        'G': set(['C', 'U']),
        'U': set(['A', 'G']),
    }

    loop_types = ['HL', 'IL']

    def current_id(self):
        with self.session() as session:
            query = session.query(mod.LoopReleases).\
                order_by(mod.LoopReleases.date.desc()).\
                limit(1)
            return query.one().loop_releases_id

    def query(self, session, pdb):
        release_id = self.current_id()
        return session.query(self.table).\
            join(mod.LoopInfo, mod.LoopInfo.loop_id == mod.LoopQa.loop_id).\
            filter(mod.LoopInfo.pdb_id == pdb).\
            filter(mod.LoopQa.release_id == release_id)

    def loops(self, pdb):
        with self.session() as session:
            lp = mod.LoopPositions
            query = session.query(mod.LoopInfo.loop_id.label('id'),
                                  mod.LoopInfo.type,
                                  mod.LoopInfo.seq,
                                  mod.UnitInfo.unit_id,
                                  mod.UnitInfo.chain).\
                join(lp, lp.loop_id == mod.LoopInfo.loop_id).\
                join(mod.UnitInfo, mod.UnitInfo == lp.unit_id).\
                filter(mod.LoopInfo.pdb_id == pdb).\
                filter(mod.LoopInfo.type.in_(self.loop_types)).\
                order_by(mod.LoopInfo.loop_id)

            loops = []
            results = it.imap(ut.row2dict, query)
            for loop_id, group in it.groupby(results, lambda r: r['id']):
                positions = list(group)
                parts = []
                loops.append({
                    'id': loop_id,
                    'sequence': positions[0]['sequence'].replace('*', ','),
                    'parts': parts,
                    'type': positions[0]['type'],
                    'chains': set(p['chain'] for p in positions),
                    'unit_ids': [p['unit_id'] for p in positions]
                })
            return loops

    def is_complementary(self, loop):
        if loop['type'] != 'IL':
            return False

        seq = loop['sequence'].replace(',', '')
        for index, char in enumerate(seq):
            rev = -1 * (index + 1)
            # print(index, char, self.complementary[char], rev, seq[rev])
            if seq[rev] not in self.complementary[char]:
                return False
        return loop['sequence']

    def has_incomplete_units(self, loop):
        with self.session() as session:
            query = session.query(mod.UnitInfo).\
                filter(mod.UnitInfo.unit_id.in_(loop['unit_ids'])).\
                filter(mod.UnitInfo.is_complete == False).\
                limit(1)
            return bool(query.count())

    def too_many_chains(self, loop):
        if loop['type'] == 'HL':
            return len(loop['chains']) == 1
        return len(loop['chains']) == 2

    def has_missing(self, loop):
        pass

    def modifications(self, loop):
        with self.session() as session:
            pos = mod.ExpSeqPosition
            umap = mod.ExpSeqUnitMapping
            query = session.query(pos.unit).\
                join(umap,
                     umap.exp_seq_position_id == pos.exp_seq_position_id).\
                filter(~pos.unit.in_(['A', 'C', 'G', 'U']))

            return [r.unit for r in query]

    def quality(self, release_id, loop):
        status = 1
        nts = ', '.join(self.get_nts(loop))
        mods = None
        compl = None

        if self.has_missing(loop):
            status = 2

        mods = self.modifications(loop)
        if mods:
            mods = ', '.join(mods)
            status = 3

        if self.too_many_chains(loop):
            status = 4

        if self.has_incomplete_units(loop):
            status = 5

        if self.is_complementary(loop):
            status = 6

        return {
            'loop_id': loop['id'],
            'status': status,
            'modifications': mods,
            'nt_signature': nts,
            'complementary': compl,
            'loop_release_id': release_id
        }

    def data(self, pdb, **kwargs):
        raise core.Skip("NYI")
        release_id = self.current_id()
        return [self.quality(release_id, loop) for loop in self.loops(pdb)]
