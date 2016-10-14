"""
"""

import operator as op
import collections as coll

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import row2dict

from pymotifs.nr import representatives as reps


class Reporter(core.Reporter):
    allow_no_data = True
    headers = [
        'Handle',
        'Version',
        'Release',
        'Release Index',
        'Method',
        'Representative',
        'Change Count'
    ]

    def release_index(self, release_id):
        with self.session() as session:
            release = session.query(mod.NrReleases).get(release_id)
            if not release:
                raise core.InvalidState("Unknown release %s" % release_id)
            return release.index

    def as_release_name(self, nr_release):
        return 'v' + nr_release

    def release_order(self, start, stop):
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
            version = rep['version']
            key = (handle, version, method)
            grouped[key][rel_name] = rep['id']

    def finalize(self, ordering, grouped):
        data = []
        for (handle, version, method), releases in grouped.items():
            last_entry = None
            entries = []
            for index, release in enumerate(ordering):
                count = 0
                if release not in releases:
                    self.logger.error("Missing release %s for %s",
                                      release, handle)
                    entries = []
                    break
                curr = releases[release]
                if index:
                    count = last_entry['Change Count']
                    prev = releases[ordering[index - 1]]
                    if prev != curr:
                        count += 1

                last_entry = {
                    'Handle': handle,
                    'Version': version,
                    'Method': method,
                    'Release': release,
                    'Release Index': index,
                    'Representative': curr,
                    'Change Count': count,
                }
                entries.append(last_entry)
            data.extend(entries)
        key = op.itemgetter('Release Index', 'Handle', 'Method')
        return sorted(data, key=key)

    def data(self, entry, **kwargs):
        start, stop, resolution = entry
        grouped = coll.defaultdict(dict)
        ordering = self.release_order(start, stop)
        for class_id in self.common_classes(start, stop, resolution):
            self.get_representatives(grouped, class_id)
        return self.finalize(ordering, grouped)
