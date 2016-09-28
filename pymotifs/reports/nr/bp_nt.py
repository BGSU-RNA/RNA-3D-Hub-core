"""This module contains logic to build a report of BP vs NT for an NR set
"""

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import row2dict


class Reporter(core.Reporter):
    headers = [
        'Group',
        'Release',
        'IFE',
        'BP',
        'NT',
        'Representative',
    ]

    def data(self, entry, **kwargs):
        release, resolution = entry
        with self.session() as session:
            query = session.query(mod.NrClasses.name.label('Group'),
                                  mod.NrClasses.nr_release_id.label('Release'),
                                  mod.NrChains.ife_id.label('IFE'),
                                  mod.IfeInfo.bp_count.label('BP'),
                                  mod.IfeInfo.length.label('NT'),
                                  mod.NrChains.rep.label('Representative'),
                                  ).\
                join(mod.NrChains,
                     mod.NrChains.nr_class_id == mod.NrClasses.nr_class_id).\
                join(mod.IfeInfo,
                     mod.IfeInfo.ife_id == mod.NrChains.ife_id).\
                filter(mod.NrClasses.nr_release_id == release).\
                filter(mod.NrClasses.resolution == resolution)
            return [row2dict(r) for r in query]
