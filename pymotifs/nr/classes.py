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

    def filter_by_resolution(self, groups, cutoff):
        if cutoff == 'all':
            return groups

        cutoff = float(cutoff)
        return groups

    def load_release(self, release_id, cutoff='4.0'):

        with self.session() as session:
            query = session.query(NrClasses.handle.label('handle'),
                                  NrClasses.version.label('version'),
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

            if cutoff != 'all':
                query = query.filter(PdbInfo.resolution <= float(cutoff))

            key = lambda v: (v.handle, v.version)
            ordered = sorted(query, key=key)
            grouped = it.groupby(ordered, key)

            results = []
            for (handle, version), members in grouped:
                results.append({
                    'members': [{'id': member.id} for member in members],
                    'handle': handle,
                    'version': int(version)
                })

        return results

    def data(self, pdbs, **kwargs):
        grouper = Grouper(self.config, self.session.maker)
        groups = grouper(pdbs, **kwargs)
        namer = Namer(self.config, self.session.maker)

        helper = Release(self.config, self.session.maker)
        release_id = helper.current('nr')

        data = []
        for named in namer(groups):
            for resolution in self.resolution_groups:
                groups = self.filter_by_resolution(named, resolution)
                naming = (resolution, named['handle'], named['version'])
                data.append(NrClasses(name='NR_%s_%s.%i' % naming,
                                      nr_release_id=release_id,
                                      resolution=resolution,
                                      handle=named['handle'],
                                      version=named['version'],
                                      comment=named['comment']))
        return data

        # Determine names for each group
        #   If unchanged reuse an old name
        #   New means a new handle
        #   Modified
        #       Only a little => bump version
        #       Largely => new handle
        #   If several parents then new name
