import itertools as it
import datetime as dt

from pymotifs import core
from pymotifs.models import NrChains
from pymotifs.models import NrClasses
from pymotifs.models import AutonomousInfo
from pymotifs.models import PdbInfo

from pymotifs.nr.release import Loader as ReleaseLoader
from pymotifs.nr.groups.simplified import Grouper
from pymotifs.nr.groups.naming import Namer
from pymotifs.utils.releases import Release
from pymotifs.utils import tmp

from pymotifs.chains.info import Loader as ChainLoader
from pymotifs.interactions import Loader as InteractionLoader


class Loader(core.MassLoader):
    dependencies = set([ReleaseLoader, ChainLoader, InteractionLoader])
    update_gap = dt.timedelta(7)  # Only update every 7 days
    resolution_groups = ['1.5', '2.0', '2.5', '3.0', '3.5', '4.0', '20.0',
                         'all']

    def has_data(self, pdbs, **kwargs):
        helper = Release(self.config, self.session.maker)
        release_id = helper.current('nr')

        with self.session() as session:
            query = session.query(NrClasses).\
                filter_by(nr_release_id=release_id)

            return bool(query.count())

    def load_release(self, release_id, cutoff='all'):
        with self.session() as session:
            query = session.query(NrClasses.handle.label('handle'),
                                  NrClasses.version.label('version'),
                                  NrClasses.id.label('class_id'),
                                  NrClasses.name.label('full_name'),
                                  AutonomousInfo.id.label('id'),
                                  ).\
                join(NrChains,
                     NrChains.nr_class_id == NrClasses.id).\
                join(AutonomousInfo,
                     AutonomousInfo.id == NrChains.autonomous_group_id).\
                join(PdbInfo,
                     PdbInfo.id == AutonomousInfo.pdb_id).\
                filter(NrClasses.nr_release_id == release_id).\
                filter(NrClasses.resolution == cutoff)

            key = lambda v: (v.handle, v.version)
            ordered = sorted(query, key=key)
            grouped = it.groupby(ordered, key)

            results = []
            for (handle, version), members in grouped:
                members = list(members)
                class_id = members[0].class_id
                full_name = members[0].full_name
                results.append({
                    'members': [{'id': member.id} for member in members],
                    'name': {
                        'class_id': class_id,
                        'full_name': full_name,
                        'handle': handle,
                        'version': int(version)
                    }
                })

        return results

    def known_handles(self):
        with self.session() as session:
            query = session.query(NrClasses.handle).distinct()
            return set(result.handle for result in query)

    def group(self, pdbs, **kwargs):
        grouper = Grouper(self.config, self.session.maker)
        return grouper(pdbs, **kwargs)

    def name(self, groups, previous):
        namer = Namer(self.config, self.session.maker)
        handles = self.known_handles()
        return namer(groups, previous, handles)

    def within_cutoff(self, chains, cutoff):
        if cutoff == 'all':
            return sorted((dict(c) for c in chains), key=lambda d: d['rank'])

        resolution = float(cutoff)
        found = [dict(c) for c in chains if c['resolution'] <= resolution]
        found.sort(key=lambda f: f['rank'])

        for index, entry in enumerate(found):
            entry['rank'] = index
        return found

    def data(self, pdbs, **kwargs):
        helper = Release(self.config, self.session.maker)
        release_id = helper.current('nr')
        previous_id = helper.previous(release_id)
        previous = self.load_release(previous_id)
        groups = self.group(pdbs, **kwargs)
        named = self.name(groups, previous)

        data = []
        mapping = {}
        for entry in named:
            for resolution in self.resolution_groups:
                parts = (resolution, entry['handle'], entry['version'])
                name = 'NR_%s_%s.%i' % parts

                filtered = self.within_cutoff(entry, resolution)
                mapping[(name, release_id)] = filtered

                data.append(NrClasses(name=name,
                                      nr_release_id=release_id,
                                      resolution=resolution,
                                      handle=entry['handle'],
                                      version=entry['version'],
                                      comment=entry['comment']))

        tmp.store('nr', mapping)
        return data
