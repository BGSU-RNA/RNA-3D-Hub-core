import abc
import collections as coll
import pprint

#import numpy as np
#import operator as op

#from sqlalchemy import func
#from sqlalchemy.orm import aliased

from pymotifs import core
from pymotifs import models as mod
from pymotifs.constants import NR_CACHE_NAME
from pymotifs.ife.info import Loader as IfeInfoLoader
from pymotifs.nr.ordering import Loader as OrderingLoader
from pymotifs.nr.representatives.using_quality import CompScore
from pymotifs.reports.nr.cqs import Groups
from pymotifs.utils import row2dict

#from .core import Representative
#from pymotifs.ife.helpers import IfeLoader


class IfeQualityLoader(core.SimpleLoader):
    """Loader to store quality data for an input equivalence class
    in table ife_cqs.
    """

    dependencies = set([IfeInfoLoader, OrderingLoader])

    """We allow this to merge data since sometimes we want to replace.

    Blake thinks that this is unnecessary/undesirable.  However, if
    a representative changes (and the IFE length as a consequence),
    all of the CQS values for that EC will need to be updated.
    """
    merge_data = True

    """We allow for no data to be written when appropriate"""
    allow_no_data = True

    """
    #Blake's sample class contained three functions, as follows:

    #def ifes(self, pdb_id):
    #    return [] # Load all IFEs for the PDB ID.

    #def ife_compscore(self, ife):
    #    return -1 # Should actually compute something

    #def data(self, pdb_id):
    #    ifes = self.ifes(pdb_id)
    #    data = []
    #    for ife in self.ifes(pdb_id): # Get all IFEs in the PDB
    #        # Let's assume that the ifes method returns a list of
    #        # dicts with the IFE information
    #        data.append(IfeCqs(
    #            ife_id=ife['id'],
    #            compscore=self.ife_compscore(ife),
    #        ))
    #    return data
    """

    # create some empty dictionaries to hold data (testing)
    hold_ri = {}
    hold_rd = {}

    def to_process(self, pdbs, **kwargs):
        """Look up the release ID to process. Ignores the given PDBs.

        Parameters
        ----------
        pdbs : list
            Ignored.

        Returns
        -------
        latest : string
            The NR release ID to process.
        """

        self.logger.info("IQL: to_process")

        latest = None
        if kwargs.get('manual', {}).get('nr_release_id', False):
            latest = kwargs['manual']['nr_release_id']
        else:
            data = self.cached(NR_CACHE_NAME)
            if not data:
                raise core.InvalidState("No precomputed grouping to store")
            latest = data['release']

        self.logger.info("IfeQualityLoader: to_process: latest: %s" % latest)

        with self.session() as session:
            return [(latest, 'foo')]


    def query(self, session, pair):
        """Create a query to calculate all entries that will be added to
        ife_cqs for the given class_id.

        Parameters
        ----------
        session : pymotifs.core.Session
            The session to use.

        pdb_id : int
            The input PDB id (to be ignored)

        Returns
        -------
        query : Query
            The query.
        """

        compscore = self._create(CompScore)
        groups = self._create(Groups)
        self.logger.info("IQL: query")

        _, class_id = pair
        self.logger.info("IQL: query: release: %s" % _)
        self.logger.info("IQL: query: class_id: %s" % class_id)

        #self.resolutions = ['1.5', '2.0', '2.5', '3.0', '3.5', '4.0', '20.0', 'all']
        #self.resolutions = ['all'] # for testing of revisions for single-pass CQS
        self.resolutions = ['1.5'] # for faster testing

        for res in self.resolutions:
            self.logger.info("IQL: query: res: %s" % res)
            #foo = groups.data(_, 1.5)
            foo = groups.data(_, res)
            for x in foo:
                self.logger.info("IQL: query: foo: %s" % x)
                #
                # need to separate "x" into two pools of data:
                # 1) resolution-independent
                # 2) resolution-dependent
                #
                # and then work out how to get those pools into the database
                #
                # Keywords:
                ## 'Compscore'
                ## 'Fraction Unobserved'
                # 'Obs Length (UI)'
                ## 'IFE ID'
                # 'Group'
                ## 'Percent Clash'
                ## 'Rfree'
                # 'Maximum Experimental Length'
                ## 'Average RSR'
                ## 'Average RSCC'
                ## 'Clashscore'
                ## 'Release'
                ## 'Resolution'
                # 'Exp Length (CI)'

                x['id'] = x['IFE ID']

                info = compscore.member_info(x)
                self.logger.info("IQL: query: foo: %s" % info)

                #localpc = compscore.percent_clash(x)

                self.hold_ri['release_id'] = _
                self.hold_rd['release_id'] = _
                self.hold_ri['ife_id'] = x['IFE ID']
                self.hold_rd['ife_id'] = x['IFE ID']
                self.hold_ri['clashscore'] = x['Clashscore']
                self.hold_ri['average_rsr'] = x['Average RSR']
                self.hold_ri['average_rscc'] = x['Average RSCC']
                self.hold_ri['percent_clash'] = x['Percent Clash']
                self.hold_ri['rfree'] = x['Rfree']
                self.hold_rd['resolution'] = x['Resolution']
                self.hold_rd['fraction_unobserved'] = x['Fraction Unobserved']
                self.hold_rd['compscore'] = x['Compscore']

                self.logger.info("IQL: query: ri_release: %s" % self.hold_ri['release_id'])
                self.logger.info("IQL: query: rd_release: %s" % self.hold_rd['release_id'])
                self.logger.info("IQL: query: ri_ife_id: %s" % self.hold_ri['ife_id'])
                self.logger.info("IQL: query: rd_ife_id: %s" % self.hold_rd['ife_id'])
                self.logger.info("IQL: query: ri_clashscore: %s" % self.hold_ri['clashscore'])
                self.logger.info("IQL: query: ri_average_rsr: %s" % self.hold_ri['average_rsr'])
                self.logger.info("IQL: query: ri_average_rscc: %s" % self.hold_ri['average_rscc'])
                self.logger.info("IQL: query: ri_percent_clash: %s" % self.hold_ri['percent_clash'])
                self.logger.info("IQL: query: ri_rfree: %s" % self.hold_ri['rfree'])
                self.logger.info("IQL: query: rd_resolution: %s" % self.hold_rd['resolution'])
                self.logger.info("IQL: query: rd_fraction_unobserved: %s" % self.hold_rd['fraction_unobserved'])
                self.logger.info("IQL: query: rd_compscore: %s" % self.hold_rd['compscore'])

            self.logger.info("IQL: query: endfoo")
        self.logger.info("IQL: query: endres")

        pass # temporary
        return session.query(mod.IfeCqs)#.\
            #filter_by(ife_id=ife_id)


    def data(self, nr_class):
        """Present only because required by the base class."""
        pass


    #def process_nr_class(self, nr_class_id):
    #    with self.session() as session:
    #        query = session.query(
    #            mod.NrChains.rank.label('index'),
    #            mod.NrChains.ife_id.label('id'),
    #            mod.IfeInfo.pdb_id,
    #            mod.IfeInfo.length,
    #            mod.IfeChains.chain_id,
    #            mod.NrClasses.name,
    #        ).\
    #            join(mod.IfeInfo, mod.IfeInfo.ife_id == mod.NrChains.ife_id).\
    #            join(mod.IfeChains,
    #                 mod.IfeChains.ife_id == mod.IfeInfo.ife_id).\
    #            join(mod.NrClasses,
    #                 mod.NrClasses.nr_class_id == mod.NrChains.nr_class_id).\
    #            filter(mod.NrClasses.nr_class_id == nr_class_id).\
    #            filter(mod.IfeChains.index == 0)

    #        data = coll.defaultdict(list)
    #        for result in query:
    #            entry = row2dict(result)
    #            entry['rep'] = (entry['index'] == 0)
    #            nr = entry['name']
    #            data[nr].append(entry)
    #    return data.values()

    pass

