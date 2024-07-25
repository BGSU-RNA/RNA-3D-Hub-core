"""
Store a new NR release and cache and NR grouping.

This will compute the next NR release id and store an entry in the NR release
table about it. It will also use the given PDBs to create a new NR grouping.
This grouping will then be cached for usage in future stages.
"""

import datetime as dt

from sqlalchemy import desc

from pymotifs import core
from pymotifs import models as mod
from pymotifs.constants import NR_CACHE_NAME

from pymotifs.utils import releases as rel
from pymotifs.nr.builder import Builder

from pymotifs.chains.info import Loader as ChainLoader
# from pymotifs.chains.species import Loader as ChainSpeciesLoader
from pymotifs.interactions.loader import Loader as InteractionLoader
from pymotifs.ife.loader import Loader as IfeLoader
from pymotifs.correspondence.loader import Loader as CorrespondenceLoader
from pymotifs.chain_chain.loader import Loader as ChainChainLoader
from pymotifs.quality.loader import Loader as QualityLoader
from pymotifs.units.loader import Loader as UnitsLoader
from pymotifs.constants import RESOLUTION_GROUPS
from pymotifs.export.ife_discrepancy import Exporter as IFEDiscrepancyExporter


class Loader(core.MassLoader):
    dependencies = set([ChainLoader, InteractionLoader,
                        IfeLoader, CorrespondenceLoader, ChainChainLoader,
                        QualityLoader, UnitsLoader, IFEDiscrepancyExporter])
    # dependencies = set([])
    update_gap = dt.timedelta(7)
    ## I think this is not safe, the following variable is just for testing purposes.
    allow_no_data = True

    def has_data(self, *args, **kwargs):
        """This will always return True because we only want to update if the time
        difference has been large enough.
        """
        return True

    def remove(self, *args, **kwargs):
        self.logger.info("Will never automatically delete nr releases"
                         " or cached data")

    def build(self, pdbs, current_release, next_release, **kwargs):
        builder = Builder(self.config, self.session)
        self.cache(NR_CACHE_NAME, builder(pdbs, current_release, next_release, cutoffs=RESOLUTION_GROUPS, **kwargs))

    def next_id(self, current):
        return rel.next_id(current, mode=self.config['release_mode']['nrlist'])

    def current_id(self):
        with self.session() as session:
            query = session.query(mod.NrReleases.nr_release_id,
                                  mod.NrReleases.index).\
                order_by(desc(mod.NrReleases.index)).\
                limit(1)

            if query.count() == 0:
                return '0.0', 0

            current = query.one()
            return current.nr_release_id, current.index

    # let's not make separate code for DNA
    # instead, let's pass in the parent and the next release number
    # we also want to make clear what type of molecule we are running on
    # def current_dna_version(self):
    #     with self.session() as session:
    #         query = session.query(mod.NrReleases.nr_release_id,
    #                               mod.NrReleases.index).\
    #             order_by(desc(mod.NrReleases.index)).\
    #             limit(1)

    #         if query.count() == 0:
    #             return '0.0', 0

    #         current = query.one()
    #         return current.nr_release_id, current.index

    def data(self, pdbs, **kwargs):

        if len(pdbs) < 100:
            raise core.Skip("Not enough pdbs to build NR release")

        now = kwargs.get('before_date', dt.datetime.strftime(dt.datetime.now(), '%Y-%m-%d %H:%M:%S'))

        if ":" in now:
            nowstring = dt.datetime.strftime(dt.datetime.strptime(now, "%Y-%m-%d %H:%M:%S").date(), "%Y%m%d")
        else:
            nowstring = now.replace("-","")

        # check command-line arguments
        nr_molecule_parent_current = kwargs.get('nr_molecule_parent_current','')
        self.logger.info("nr_molecule_parent_current: %s" % nr_molecule_parent_current)

        if nr_molecule_parent_current:
            # manually set the molecule, parent release, current release
            fields = nr_molecule_parent_current.split(",")
            molecule = fields[0]   # DNA or RNA
            parent = fields[1]     # like 0.3
            next = fields[2]       # like 0.4
            if parent == '0.0':
                parent = next
            pass
        else:
            current, index = self.current_id()
            next = self.next_id(current)
            parent = current
            if current == '0.0':
                parent = next

        self.logger.info('Building NR release %s with parent %s' % (next, parent))

        # build the current release; data must be cached and will be read in a later stage
        self.build(pdbs, parent, next, **kwargs)

        self.logger.info('Built NR release %s, saving to nr_releases table' % next)

        if nr_molecule_parent_current:
            raise core.Skip("Filling in DNA releases, no need to write to nr_releases table")
        else:
            # regular new release
            return mod.NrReleases(nr_release_id=next,
                                date=now,
                                parent_nr_release_id=parent,
                                description = nowstring,
                                index=index + 1)

