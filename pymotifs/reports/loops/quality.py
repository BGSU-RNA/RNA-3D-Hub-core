"""
"""

import operator as op
import itertools as it
import collections as coll

from pymotifs import core
from pymotifs import models as mod

from pymotifs.utils import row2dict
from pymotifs.motifs.release import Loader


class Reporter(core.Reporter):
    headers = [
        'Loop',
        'Type',
        'Motif',
        'Nt',
        'Internal Pairs',
        'External Pairs',
        'Internal Stacks',
        'External Stacks',
        'Internal Basephosphate',
        'External Basephosphate',
        'Internal Interactions',
        'External Interactions',
        'Total Interactions',
        'RSR',
        'RSRZ',
    ]

    def group(self, results):
        grouped = it.groupby(results, op.attrgetter('pdb_id'))
        data = []
        for pdb, pairs in grouped:
            data.append((pdb, list(it.imap(op.itemgetter(1), pairs))))
        return data

    def loops_in_pdbs(self, pdbs):
        with self.session() as session:
            query = session.query(mod.LoopInfo.pdb_id,
                                  mod.LoopInfo.loop_id,
                                  ).\
                filter(mod.LoopInfo.pdb_id.in_(pdbs)).\
                distinct().\
                order_by(mod.LoopInfo.pdb_id)
            return self.group(query)

    def loops_in_motifs(self, release):
        with self.session() as session:
            query = session.query(mod.MlLoops.loop_id,
                                  mod.LoopInfo.pdb_id,
                                  ).\
                filter(mod.MlLoops.ml_release_id == release).\
                distinct().\
                order_by(mod.LoopInfo.pdb_id)
            return self.group(query)

    def loops_in_nr(self, release, resolution='all', **kwargs):
        with self.session() as session:
            ifes = mod.IfeInfo
            pos = mod.LoopPositions
            units = mod.UnitInfo
            chains = mod.ChainInfo
            ife_chains = mod.IfeChains
            nr = mod.NrChains
            classes = mod.NrClasses
            loops = mod.LoopInfo
            query = session.query(loops.loop_id,
                                  loops.pdb_id,
                                  ).\
                join(pos, pos.loop_id == loops.loop_id).\
                join(units, units.unit_id == pos.unit_id).\
                join(chains, (chains.pdb_id == units.pdb_id) &
                     (chains.chain_name == units.chain)).\
                join(ife_chains, ife_chains.chain_id == chains.chain_id).\
                join(ifes, ifes.ife_id == ife_chains.ife_id).\
                join(nr, nr.ife_id == ifes.ife_id).\
                join(classes, nr.nr_class_id == classes.nr_class_id).\
                filter(nr.nr_release_id == release).\
                filter(nr.rep == 1).\
                filter(classes.resolution == resolution).\
                distinct().\
                order_by(mod.LoopInfo.pdb_id)
            return self.group(query)

    def all_loops(self):
        with self.session() as session:
            query = session.query(mod.LoopInfo.pdb_id,
                                  mod.LoopInfo.loop_id,
                                  ).\
                distinct().\
                order_by(mod.LoopInfo.pdb_id)
            return self.group(query)

    def to_process(self, pdbs, nr_release=None, motif_release=None, **kwargs):
        if nr_release:
            loops = self.loops_in_nr(nr_release, **kwargs)
        if motif_release:
            loops = self.loops_in_motif(motif_release)
        if pdbs:
            loops = self.loops_in_pdbs(pdbs)
        else:
            loops = self.all_loops()
        return [tuple(loops)]

    def current_motif_release(self):
        loader = Loader(self.config, self.session)
        release, _ = loader.current_id()
        return release

    def pairs(self, pdb):
        pairs = coll.defaultdict(lambda: {'Pairs': set(), 'Stacks': set(),
                                          'Basephosphate': set()})
        with self.session() as session:
            interactions = mod.UnitPairsInteractions
            query = session.query(interactions.unit_id_1,
                                  interactions.unit_id_2,
                                  interactions.f_lwbp.label('Pairs'),
                                  interactions.f_stacks.label('Stacks'),
                                  interactions.f_bphs.label('Basephosphate'),
                                  ).\
                filter(interactions.pdb_id == pdb)
            for result in query:
                data = row2dict(result)
                unit1 = data.pop('unit_id_1')
                unit2 = data.pop('unit_id_2')
                if unit1 == unit2:
                    continue
                for name, value in data.items():
                    if value and not value.startswith('n'):
                        pairs[unit1][name].add(unit2)
        return pairs

    def positions(self, pdb):
        positions = coll.defaultdict(set)
        with self.session() as session:
            info = mod.LoopInfo
            pos = mod.LoopPositions
            query = session.query(pos.loop_id,
                                  pos.unit_id,
                                  ).\
                join(info, info.loop_id == pos.loop_id).\
                filter(info.pdb_id == pdb)
            for result in query:
                positions[result.loop_id].add(result.unit_id)
        return positions

    def loop_quality(self, pdb, loop_ids, motif_release=None, **kwargs):
        if not motif_release:
            release_id = self.current_motif_release()

        pairs = self.pairs(pdb)
        known_positions = self.positions(pdb)
        with self.session() as session:
            info = mod.LoopInfo
            positions = mod.LoopPositions
            quality = mod.UnitQuality
            motifs = mod.MlLoops
            query = session.query(info.loop_id.label('Loop'),
                                  info.type.label('Type'),
                                  motifs.motif_id.label('Motif'),
                                  positions.unit_id.label('Nt'),
                                  quality.real_space_r.label('RSR'),
                                  quality.z_score.label('RSRZ'),
                                  ).\
                join(positions, positions.loop_id == info.loop_id).\
                outerjoin(motifs, (motifs.loop_id == info.loop_id) &
                          (motifs.ml_release_id == release_id)).\
                outerjoin(quality, quality.unit_id == positions.unit_id).\
                filter(info.pdb_id == pdb).\
                filter(info.loop_id.in_(loop_ids))

            for result in query:
                entry = row2dict(result)
                nts = known_positions[entry['Loop']]
                entry['Internal Pairs'] = 0
                entry['External Pairs'] = 0
                entry['Internal Stacks'] = 0
                entry['External Stacks'] = 0
                entry['Internal Basephosphate'] = 0
                entry['External Basephosphate'] = 0
                entry['Total Interactions'] = 0
                for nt in nts:
                    for key, curr in pairs[nt].items():
                        entry['Internal ' + key] = len(curr.intersection(nts))
                        entry['External ' + key] = len(curr - nts)
                entry['Internal Interactions'] = entry['Internal Pairs'] + \
                    entry['Internal Stacks'] + entry['Internal Basephosphate']
                entry['External Interactions'] = entry['External Pairs'] + \
                    entry['External Stacks'] + entry['External Basephosphate']
                entry['Total Interactions'] += entry['Internal Interactions']
                entry['Total Interactions'] += entry['External Interactions']
                yield entry

    def data(self, loops, **kwargs):
        data = []
        for (pdb, loops) in loops:
            data.extend(self.loop_quality(pdb, loops, **kwargs))
        return data
