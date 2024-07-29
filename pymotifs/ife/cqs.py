"""This module contains the logic to create IFE-level CQS
(Composite Quality Scoring) data for subsequent import
into the database.
"""

import operator as op

from sqlalchemy import func
from sqlalchemy.orm import aliased

from copy import deepcopy

from pymotifs import core
from pymotifs import models as mod
from pymotifs.ife.helpers import IfeLoader
from pymotifs.ife.info import Loader as IfeInfoLoader
from pymotifs.quality.loader import Loader as QualityLoader
from pymotifs.exp_seq.loader import Loader as ExpSeqLoader
from pymotifs.nr.representatives.using_quality import CompScore
from pymotifs.utils import row2dict


class IfeQualityLoader(core.SimpleLoader):
    """Loader to store non-release-dependent (i.e., IFE-based)
    quality data for an input equivalence class in table
    ife_cqs.
    """

    dependencies = set([IfeInfoLoader,QualityLoader,ExpSeqLoader])

    """Handle replacements via the recalculate option."""
    merge_data = True

    """We allow for no data to be written when appropriate."""
    allow_no_data = True

    mark = False  # do not try to mark in pdb_analysis_status

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

    def to_process(self, pdbs, **kwargs):
        """
        Look up the list of IFEs to process, from among pdb ids in pdbs

        Parameters
        ----------
        pdbs : list
            Ignored.

        Returns
        -------
        release : string
            The NR release ID to process.
        """

        resolution = 'all'
        release = None

        pdbs_set = set(pdbs)
        ife_list = []

        if False and kwargs.get('manual', {}).get('nr_release_id', False):
            # This query would not need to rewritten to get the release id
            # from nr_classes and join with this nr_chains table
            release = kwargs['manual']['nr_release_id']
            self.logger.info("IQL: to_process: release: %s" % release)
            with self.session() as session:
                query = session.query(mod.NrChains.ife_id).\
                    join(mod.NrClasses,
                         mod.NrClasses.nr_class_id == mod.NrChains.nr_class_id).\
                    filter(mod.NrClasses.nr_release_id == release).\
                    filter(mod.NrClasses.resolution == resolution)
                return [r.ife_id for r in query]
        else:
            with self.session() as session:
                query = session.query(mod.IfeInfo.ife_id).\
                    filter(mod.IfeInfo.model.isnot(None)).\
                    filter(mod.IfeInfo.new_style == True)
                for r in query:
                    pdb_id = r.ife_id.split('|')[0]
                    if pdb_id in pdbs_set:
                        ife_list.append(r.ife_id)

                if len(ife_list) == 0:
                    raise core.Skip("No IFEs to process for CQS")

                return ife_list

        pass

    def query(self, session, ife_id):
        return session.query(mod.IfeCqs.ife_id).filter_by(ife_id=ife_id)


    def class_property(self, ifes, name):
        return {ife[name] for ife in ifes}

    def ife_info(self, nr_class):
        ife_id = self.class_property(nr_class, 'id')
        with self.session() as session:
            data = {}

            query = session.query(mod.IfeInfo.ife_id).\
                filter(mod.IfeInfo.ife_id == ife_id)

            for result in query:
                entry = row2dict(result)
                ife_id = entry.pop('ife_id')
                chain_ids = ife_id.split('+')
                chains = [p.split('|')[-1] for p in chain_ids]
                entry['Chains'] = ', '.join(chains)
                data[ife_id] = entry
        return data

    # comment this out on 10/13/2023
    # we cannot find anywhere that uses the following functions.
    # If these two functions need to be used in the future, make the query works corrently with NrChains table/release_id.
    # def nr_classes(self, release, resolution):
    #     return self.load_nr_classes(release, resolution)

    # def load_nr_classes(self, release, resolution):
    #     with self.session() as session:
    #         query = session.query(
    #             mod.NrChains.rank.label('index'),
    #             mod.NrChains.ife_id.label('id'),
    #             mod.IfeInfo.pdb_id,
    #             mod.IfeInfo.length,
    #             mod.IfeChains.chain_id,
    #             mod.NrClasses.name,
    #         ).\
    #             join(mod.IfeInfo, mod.IfeInfo.ife_id == mod.NrChains.ife_id).\
    #             join(mod.IfeChains,
    #                  mod.IfeChains.ife_id == mod.IfeInfo.ife_id).\
    #             join(mod.NrClasses,
    #                  mod.NrClasses.nr_class_id == mod.NrChains.nr_class_id).\
    #             filter(mod.NrClasses.nr_release_id == release).\
    #             filter(mod.NrClasses.resolution == resolution).\
    #             filter(mod.IfeChains.index == 0)

    #         data = coll.defaultdict(list)
    #         for result in query:
    #             entry = row2dict(result)
    #             entry['rep'] = (entry['index'] == 0)
    #             nr = entry['name']
    #             data[nr].append(entry)
    #     return data.values()

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
                'Resolution': quality['resolution'],
            }
        return data

    def member_info(self, member):
        with self.session() as session:
            info = session.query(mod.IfeInfo.pdb_id.label('pdb'),
                                 mod.IfeInfo.model).\
                filter_by(ife_id=member['id']).\
                one()
            info = row2dict(info)
            info.update(member)

            with self.session() as session:
                query = session.query(mod.ChainInfo.chain_name,
                                      mod.IfeChains.is_structured,
                                      ).\
                    join(mod.IfeChains,
                         mod.IfeChains.chain_id == mod.ChainInfo.chain_id).\
                    filter_by(ife_id=member['id'])

                if not query.count():
                    raise core.InvalidState("Could not find chains for %s" %
                                            member)

                all_chains = [row2dict(c) for c in query]
                valid = op.itemgetter('is_structured')
                chains = [c['chain_name'] for c in all_chains if valid(c)]
                if not chains:
                    chains = [c['chain_name'] for c in all_chains]

            info['chains'] = chains
            loader = self._create(IfeLoader)
            info['sym_op'] = loader.sym_op(info['pdb'])

            return info

    def observed_length(self, info):
        self.logger.info("observed_length: info: %s" % info)
        with self.session() as session:
            query = session.query(mod.UnitInfo.chain_index).\
                distinct().\
                filter(mod.UnitInfo.pdb_id == info['pdb']).\
                filter(mod.UnitInfo.model == info['model']).\
                filter(mod.UnitInfo.sym_op == info['sym_op']).\
                filter(mod.UnitInfo.chain.in_(info['chains'])).\
                filter(mod.UnitInfo.unit.in_(['A', 'C', 'G', 'U'])).\
                filter(mod.UnitInfo.chain_index != None)
            self.logger.info("observed_length: after query generation")
            return query.count()


    def data(self, entry, **kwargs):
        """Create a report about the NR set.

        Parameters
        ----------
        entry : str
            The IFE for which to collect IFE-level composite
            quality score data.

        Returns
        -------
            The required data for the database update step.
        """

        ife = dict([('index', 0), ('id', entry)])
        ife = self.member_info(ife)
        ife['length'] = self.observed_length(ife)

        nr_class = []
        nr_class.append(ife)

        ife_info = self.ife_info(nr_class)
        quality_data = self.ife_quality_data(nr_class)

        data = {
            'IFE ID': ife['id'],
            'Observed Length': ife['length'],
        }
        data.update(ife_info[ife['id']])
        data.update(quality_data[ife['id']])
        self.logger.debug("data: data: %s" % data)

        yield mod.IfeCqs(
            ife_id = data['IFE ID'],
            obs_length = data['Observed Length'],
            clashscore = data['Clashscore'],
            average_rsr = data['Average RSR'],
            average_rscc = data['Average RSCC'],
            percent_clash = data['Percent Clash'],
            rfree = data['Rfree'],
            resolution = data['Resolution'])

