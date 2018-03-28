"""This module contains the logic to create IFE-level CQS
(Composite Quality Scoring) data for subsequent import 
into the database.
"""

import abc
import collections as coll
import pprint

#import itertools as it
#import numpy as np
#import operator as op

#from sqlalchemy import func
#from sqlalchemy.orm import aliased
#from sqlalchemy.sql import operators
#from sqlalchemy.sql.expression import func
#from sqlalchemy.sql.functions import coalesce

from copy import deepcopy

from pymotifs import core
from pymotifs import models as mod
from pymotifs.constants import NR_CACHE_NAME
from pymotifs.ife.info import Loader as IfeInfoLoader
#from pymotifs.nr.ordering import Loader as OrderingLoader
from pymotifs.nr.representatives.using_quality import CompScore
from pymotifs.utils import row2dict

#from .core import Representative
#from pymotifs.ife.helpers import IfeLoader

#operators._PRECEDENCE['SEPARATOR'] = 0


class IfeQualityLoader(core.SimpleLoader):
    """Loader to store non-release-dependent (i.e., IFE-based) 
    quality data for an input equivalence class in table
    ife_cqs.
    """

    dependencies = set([IfeInfoLoader])
    #dependencies = set([IfeInfoLoader, OrderingLoader])

    """Handle replacements via the recalculate option."""
    merge_data = True

    """We allow for no data to be written when appropriate."""
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

    def data(self, entry, **kwargs):
        """Create a report about the NR set.

        Parameters
        ----------
        release : str
            The NR release to report on.

        resolution : str
            The resolution to use.

        entry : tuple
            ( release, resolution )

        Returns
        -------
        rows : list
            A list of dicts to write for the report.
        """

        release, resolution = entry
        nr_classes = self.nr_classes(release, resolution)
        index = 0
        for nr_class in nr_classes:
            ife_info = self.ife_info(nr_class)
            quality_data = self.ife_quality_data(nr_class)
            for ife in nr_class:
                self.logger.debug("Processing ife: %s" % ife['id'])
                data = {
                    'Release': release,
                    'IFE ID': ife['id'],
                }
                data.update(ife_info[ife['id']])
                data.update(quality_data[ife['id']])
                yield data
                index + 1

    def class_property(self, ifes, name):
        return {ife[name] for ife in ifes}

    def ife_info(self, ifes):
        ife_ids = self.class_property(ifes, 'id')
        with self.session() as session:
            query = session.query(mod.IfeInfo.ife_id).\
                filter(mod.IfeInfo.ife_id.in_(ife_ids))

            data = {}
            for result in query:
                entry = row2dict(result)
                ife_id = entry.pop('ife_id')
                chain_ids = ife_id.split('+')
                chains = [p.split('|')[-1] for p in chain_ids]
                entry['Chains'] = ', '.join(chains)
                data[ife_id] = entry
        return data

    def nr_classes(self, release, resolution):
        return self.load_nr_classes(release, resolution)

    def load_nr_classes(self, release, resolution):
        with self.session() as session:
            query = session.query(
                mod.NrChains.rank.label('index'),
                mod.NrChains.ife_id.label('id'),
                mod.IfeInfo.pdb_id,
                mod.IfeInfo.length,
                mod.IfeChains.chain_id,
                mod.NrClasses.name,
            ).\
                join(mod.IfeInfo, mod.IfeInfo.ife_id == mod.NrChains.ife_id).\
                join(mod.IfeChains,
                     mod.IfeChains.ife_id == mod.IfeInfo.ife_id).\
                join(mod.NrClasses,
                     mod.NrClasses.nr_class_id == mod.NrChains.nr_class_id).\
                filter(mod.NrClasses.nr_release_id == release).\
                filter(mod.NrClasses.resolution == resolution).\
                filter(mod.IfeChains.index == 0)

            data = coll.defaultdict(list)
            for result in query:
                entry = row2dict(result)
                entry['rep'] = (entry['index'] == 0)
                nr = entry['name']
                data[nr].append(entry)
        return data.values()

    def ife_quality_data(self, ifes):
        compscore = self._create(CompScore)
        members = deepcopy(ifes)
        compscore.load_quality(members)
        data = {}
        for ife in members:
            quality = ife['quality']
            self.logger.debug('quality: %s' % quality)
            data[ife['id']] = {
                'Clashscore': quality.get('clashscore', 100),
                'Average RSR': quality['average_rsr'],
                'Percent Clash': quality['percent_clash'],
                'Average RSCC': quality['average_rscc'],
                'Rfree': quality['rfree'],
            }
        return data

    def to_process(self, pdbs, **kwargs):
        """Look up the release ID to process, which the next step will ignore.
        Ignores the given PDBs.

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

        self.logger.info("IQL: to_process: latest: %s" % latest)

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

        self.logger.info("IQL: query")

        release, class_id = pair
        resolution = 'all'
        resolution = '1.5' # for faster testing -- delete when finished

        self.logger.info("IQL: query: release: %s" % release)
        self.logger.info("IQL: query: resolution: %s" % resolution)
        self.logger.info("IQL: query: class_id: %s" % class_id)

        cqs_data = {} 
        self.logger.info("IQL: query: cqs_data (in): %s" % cqs_data) 
        cqs_data = self.data((release, resolution)) 
        for item in cqs_data:
            self.logger.info("IQL: query: cqs_data (out): %s" % item) 

        pass # temporary
        return session.query(mod.IfeCqs)#.\
            #filter_by(ife_id=ife_id)

    pass

