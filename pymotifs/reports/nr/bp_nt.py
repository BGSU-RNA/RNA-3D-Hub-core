"""This module contains logic to build a report of BP vs NT for an NR set
"""

import operator as op
import itertools as it
import collections as coll

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import row2dict

from pymotifs.nr import representatives as reps


class Reporter(core.Reporter):

    @property
    def headers(self):
        known = sorted(name for name, _ in reps.known())
        return [
            'Group',
            'Release',
            'IFE',
            'BP',
            'NT',
            'Ratio',
            'Current',
        ] + known

    def chains(self, release_id, resolution):
        with self.session() as session:
            chains = mod.NrChains
            ife = mod.IfeInfo
            pdbs = mod.PdbInfo
            classes = mod.NrClasses
            query = session.query(chains.nr_release_id,
                                  classes.name,
                                  classes.handle,
                                  classes.version,
                                  ife.ife_id.label('id'),
                                  ife.bp_count.label('bp'),
                                  ife.length,
                                  pdbs.resolution,
                                  pdbs.experimental_technique.label('method'),
                                  ).\
                join(ife, ife.ife_id == chains.ife_id).\
                join(classes, classes.nr_class_id == chains.nr_class_id).\
                join(pdbs, pdbs.pdb_id == ife.pdb_id).\
                filter(classes.nr_release_id == release_id).\
                filter(classes.resolution == resolution).\
                order_by(classes.name)

            found = [row2dict(r) for r in query]
            grouped = it.groupby(found, op.itemgetter('name'))
            return [list(g) for n, g in grouped]

    def representatives(self, release_id, resolution, **kwargs):
        def empty():
            known = [name for name, _ in reps.known()]
            return dict(zip(known, [False] * len(known)))

        data = coll.defaultdict(empty)
        for chains in self.chains(release_id, resolution):
            for method, klass in reps.known():
                finder = klass(self.config, self.session)
                rep = finder(chains)
                data[rep['id']][method] = True
        return data

    def chain_status(self, reps, release, resolution, **kwargs):
        with self.session() as session:
            query = session.query(mod.NrClasses.name.label('Group'),
                                  mod.NrClasses.nr_release_id.label('Release'),
                                  mod.NrChains.ife_id.label('IFE'),
                                  mod.IfeInfo.bp_count.label('BP'),
                                  mod.IfeInfo.length.label('NT'),
                                  mod.NrChains.rep,
                                  ).\
                join(mod.NrChains,
                     mod.NrChains.nr_class_id == mod.NrClasses.nr_class_id).\
                join(mod.IfeInfo,
                     mod.IfeInfo.ife_id == mod.NrChains.ife_id).\
                filter(mod.NrClasses.nr_release_id == release).\
                filter(mod.NrClasses.resolution == resolution)

            data = []
            for result in query:
                entry = row2dict(result)
                rep = entry.pop('rep')
                entry['Current'] = 'Member'
                entry['Ratio'] = round(float(entry['BP']) / entry['NT'], 4)
                if rep:
                    entry['Current'] = 'Representative'
                for method, status in reps[entry['IFE']].items():
                    entry[method] = 'Member'
                    if status:
                        entry[method] = 'Representative'
                data.append(entry)
            return data

    def data(self, entry, **kwargs):
        release, resolution = entry
        reps = self.representatives(release, resolution)
        return self.chain_status(reps, release, resolution, **kwargs)
