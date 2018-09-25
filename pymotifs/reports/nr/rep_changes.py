"""This module contains the logic to build a report to describe the
"""

import operator as op
import collections as coll

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import row2dict

from pymotifs.nr import representatives as reps


class Reporter(core.Reporter):
    allow_no_data = True

    def create_headers(self, start, stop, resolution):
        headers = ['Handle', 'Updates', 'Method']
        headers.extend(self.release_names(start, stop, resolution))
        headers.append('Changes')
        return headers

    def release_index(self, release_id):
        with self.session() as session:
            release = session.query(mod.NrReleases).get(release_id)
            if not release:
                raise core.InvalidState("Unknown release %s" % release_id)
            return release.index

    def as_release_name(self, nr_release):
        return 'v' + nr_release

    def release_names(self, start, stop, resolution):
        start_index = self.release_index(start)
        stop_index = self.release_index(stop)
        with self.session() as session:
            query = session.query(mod.NrReleases.nr_release_id).\
                filter(mod.NrReleases.index >= start_index).\
                filter(mod.NrReleases.index <= stop_index).\
                order_by(mod.NrReleases.index)
            return [self.as_release_name(r.nr_release_id) for r in query]

    def possible_classes(self, start, stop, resolution):
        start_index = self.release_index(start)
        stop_index = self.release_index(stop)
        rel = mod.NrReleases
        with self.session() as session:
            query = session.query(mod.NrClasses.handle,
                                  mod.NrClasses.nr_release_id,
                                  mod.NrClasses.nr_class_id,
                                  ).\
                join(rel,
                     mod.NrClasses.nr_release_id == rel.nr_release_id).\
                filter(rel.index >= start_index).\
                filter(rel.index <= stop_index).\
                filter(mod.NrClasses.resolution == resolution).\
                order_by(mod.NrClasses.nr_release_id)
            return [row2dict(r) for r in query]

    def common_classes(self, start, stop, resolution):
        possible = self.possible_classes(start, stop, resolution)
        releases = {p['nr_release_id'] for p in possible}
        found = coll.defaultdict(set)
        for entry in possible:
            found[entry['handle']].add(entry['nr_release_id'])

        common = set()
        for handle, present_in in found.items():
            if present_in == releases:
                common.add(handle)

        return [e['nr_class_id'] for e in possible if e['handle'] in common]

    def chains(self, class_id):
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
                filter(classes.nr_class_id == class_id)

            return [row2dict(r) for r in query]

    def methods(self, methods='*', **kwargs):
        if methods == '*':
            return reps.known()
        selected = []
        methods = set(methods)
        for name, klass in reps.known():
            if name in methods:
                selected.append((name, klass))
        return selected

    def get_representatives(self, grouped, class_id, **kwargs):
        chains = self.chains(class_id)
        for method, klass in self.methods(**kwargs):
            finder = klass(self.config, self.session)
            rep = finder(chains)
            rel_name = self.as_release_name(rep['nr_release_id'])
            handle = rep['handle']
            version = int(rep['version'])
            key = (handle, method)
            grouped[key][rel_name] = rep['id']
            grouped[key]['Updates'].append(version)

    def finalize(self, grouped):
        data = []
        for (handle, method), releases in grouped.items():
            versions = releases.pop('Updates')
            entry = {
                'Handle': handle,
                'Method': method,
                'Updates': max(versions) - min(versions),
            }
            entry.update(releases)
            entry['Changes'] = len(set(releases.values()))
            data.append(entry)
        return sorted(data, key=op.itemgetter('Handle', 'Method'))

    def data(self, entry, **kwargs):
        start, stop, resolution = entry
        self.headers = self.create_headers(start, stop, resolution)
        grouped = coll.defaultdict(lambda: {'Updates': []})
        for class_id in self.common_classes(start, stop, resolution):
            self.get_representatives(grouped, class_id)
        return self.finalize(grouped)
