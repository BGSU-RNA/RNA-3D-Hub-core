"""This module contains the logic to create NR-level CQS (Composite
Quality Score) data for subsequent import into the database.  Uses
IFE-specific data that were previously generated using ife.cqs and
stored in the database (table ife_cqs).
"""

import abc
import collections as coll
import pprint

from pymotifs import core
from pymotifs import models as mod
from pymotifs.constants import COMPSCORE_COEFFICENTS
from pymotifs.constants import NR_CACHE_NAME
# from pymotifs.nr.chains import Loader as ChainLoader
from pymotifs.nr.class_rank import Loader as ClassRankLoader
from pymotifs.nr.parent_counts import Loader as CountLoader
from pymotifs.utils import row2dict


class NrQualityLoader(core.SimpleLoader):
    """Loader to store quality data for an input equivalence class
    in table nr_cqs.
    """

    dependencies = set([ClassRankLoader, CountLoader])

    """We allow this to merge data since sometimes we want to replace.

    Blake thinks that this is unnecessary/undesirable.  However, if
    class membership changes (and the maximum IFE length as a
    consequence), all of the CQS values for that EC will need to be
    updated.
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

    def to_process(self, pdbs, **kwargs):
        """Collect the list of nr_class name values to process for the
        specified release. Ignores the given PDBs.

        Parameters
        ----------
        pdbs : list
            Ignored.

        Returns
        -------
        classlist : list
            The list of NR class names to process.
        """

        resolution = 'all'

        latest = None
        if kwargs.get('manual', {}).get('nr_release_id', False):
            latest = kwargs['manual']['nr_release_id']
        else:
            data = self.cached(NR_CACHE_NAME)
            if not data:
                raise core.InvalidState("No precomputed grouping to store")
            latest = data['release']

        classlist = self.list_nr_classes(latest, resolution)

        with self.session() as session:
            return classlist


    def load_ife_cqs_data(self, ife_list, nr_name):
        with self.session() as session:
            query = session.query(
                mod.IfeCqs.ife_id,
                mod.IfeCqs.obs_length,
                mod.IfeCqs.clashscore,
                mod.IfeCqs.average_rsr,
                mod.IfeCqs.average_rscc,
                mod.IfeCqs.percent_clash,
                mod.IfeCqs.rfree,
                mod.IfeCqs.resolution,
                ).\
                filter(mod.IfeCqs.ife_id.in_(ife_list))

            data = coll.defaultdict(list)

            max_exp_len = 0

            for result in query: 
                entry = row2dict(result)
                ii = entry['ife_id']
                entry['nr_name'] = nr_name
                data[ii].append(entry)
                if result[1] > max_exp_len:
                    max_exp_len = result[1]

        for ife in ife_list:
            if data[ife]:
                ife_data = data[ife]
                obs_length = ife_data[0]['obs_length']
                ife_data[0]['max_exp_len'] = max_exp_len
            else:
                self.logger.warning("NQL: data: LICD: no data for %s" % ife)
                continue
            truth, fraction_unobserved = self.fraction_unobserved(obs_length, max_exp_len)
            percent_observed = (1 - fraction_unobserved)
            data[ife][0]['fraction_unobserved'] = fraction_unobserved
            data[ife][0]['percent_observed'] = percent_observed
            compscore = self.compscore(data[ife])
            data[ife][0]['compscore'] = compscore

        return data.values()

    def list_nr_classes(self, release, resolution):
        with self.session() as session:
            query = session.query(mod.NrClasses.name).\
                filter(mod.NrClasses.nr_release_id == release).\
                filter(mod.NrClasses.resolution == resolution)

            data = []
            for result in query:
                data.append(result[0])
        return data

    def query(self, session, nr_name):
        return session.query(mod.NrCqs.nr_name).\
            filter(mod.NrCqs.nr_name == nr_name)

    def data(self, nr_name, **kwargs):
        """Collect composite quality scoring data for the NR class.

        Parameters
        ----------
        nr_name : str
            Name of class for which to collect NR-level composite
            quality score (CQS) data.

        Returns
        -------
            The required data for the database update step.
        """

        cqs_data = {}

        ife_list = []
        # we stop using the following query cuz we stop updating the nrchains table
        # with self.session() as session:
        #     query = session.query(mod.NrChains.ife_id).\
        #         join(mod.NrClasses, mod.NrChains.nr_class_id == mod.NrClasses.nr_class_id).\
        #         filter(mod.NrClasses.name == nr_name)
        with self.session() as session:
            query = session.query(mod.NrClassRank.ife_id).\
                filter(mod.NrClassRank.nr_class_name == nr_name)

            for result in query:
                ife_list.append(result[0])
    
        data = self.load_ife_cqs_data(ife_list, nr_name)

        for ife_output in data:
            for ife_out in ife_output:
                yield mod.NrCqs(
                    ife_id = ife_out['ife_id'],
                    nr_name = nr_name,
                    maximum_experimental_length = ife_out['max_exp_len'],
                    fraction_unobserved = ife_out['fraction_unobserved'],
                    percent_observed = ife_out['percent_observed'],
                    composite_quality_score = ife_out['compscore'])

    def fraction_unobserved(self, obs, max):
        observed = float(obs)
        experimental = float(max)
        if experimental == 0:
            return (True, 1)
        return (True, (1 - (observed / experimental)))

    def compscore(self, member):
        """
        Compute composite quality score using six indicators weighted by
        various coefficients set in constants.py
        """
        compscore = COMPSCORE_COEFFICENTS['resolution'] * member[0]['resolution']
        compscore += COMPSCORE_COEFFICENTS['percent_clash'] * member[0]['percent_clash']
        compscore += COMPSCORE_COEFFICENTS['average_rsr'] * member[0]['average_rsr']
        compscore += COMPSCORE_COEFFICENTS['average_rscc'] * (1 - member[0]['average_rscc'])
        compscore += COMPSCORE_COEFFICENTS['rfree'] * member[0]['rfree']
        compscore += COMPSCORE_COEFFICENTS['fraction_unobserved'] * member[0]['fraction_unobserved']

        if compscore < 0:
            raise core.InvalidState("Invalid compscore (%s) for %s" % (compscore, member))

        return compscore

    pass

