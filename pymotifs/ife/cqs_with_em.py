"""
This stage calculates the components that are used in nr_cqs to
compute the Composite Quality Score for IFEs
"""

import numpy as np
import operator as op

from sqlalchemy import func
from sqlalchemy.orm import aliased


from pymotifs import core
from pymotifs import models as mod
from pymotifs.ife.info import Loader as IfeInfoLoader
from pymotifs.quality.loader import Loader as QualityLoader
from pymotifs.exp_seq.loader import Loader as ExpSeqLoader
from pymotifs.utils import row2dict


class IfeQualityLoader(core.SimpleLoader):
    """
    Loader to store non-release-dependent (i.e., IFE-based)
    quality data for an input equivalence class in table
    ife_cqs.
    CQS also requires fraction observed, which changes
    according to the other structures in the equivalence class.
    """

    dependencies = set([IfeInfoLoader,QualityLoader,ExpSeqLoader])

    """Handle replacements via the recalculate option."""
    merge_data = True

    """We allow for no data to be written when appropriate."""
    allow_no_data = True

    mark = False  # do not try to mark in pdb_analysis_status

    testing = False
    fixing = False

    def to_process(self, pdbs, **kwargs):
        """
        Look up the list of IFEs to process, from among pdb ids in pdbs

        Parameters
        ----------
        pdbs : list

        Returns
        -------
        ife_list : list
            List of IFEs to process
        """

        # get list of all current ifes in the given pdbs
        with self.session() as session:
            query = session.query(mod.IfeInfo.ife_id).\
                filter(mod.IfeInfo.model.isnot(None)).\
                filter(mod.IfeInfo.new_style == True).\
                filter(mod.IfeInfo.pdb_id.in_(pdbs))

            ife_list = set([r.ife_id for r in query])

        if len(pdbs) == 1:
            return ife_list

        if len(ife_list) == 0:
            raise core.Skip("No IFEs to process for CQS")

        # get list of ifes that already have entries in ife_cqs table
        with self.session() as session:
            query = session.query(mod.IfeCqs.ife_id)

            ife_processed = set([r.ife_id for r in query])

        # identify ifes that were never processed
        ife_never_processed = ife_list - ife_processed

        # some EM structures have not yet had average_Q_score computed

        # look at pdbs solved by EM
        with self.session() as session:
            query = session.query(mod.PdbInfo.pdb_id).\
                filter(mod.PdbInfo.experimental_technique == "ELECTRON MICROSCOPY").\
                filter(mod.PdbInfo.pdb_id.in_(pdbs))
            em_pdbs = set([row.pdb_id for row in query])
        self.logger.info('Found %d EM pdbs' % len(pdbs))

        # find ifes that have average_Q_score
        with self.session() as session:
            query = session.query(mod.IfeCqs.ife_id).\
                        filter(mod.IfeCqs.average_Q_score.isnot(None))

            ife_has_Q_score = set([r.ife_id for r in query])

        ife_no_Q_score = ife_list - ife_has_Q_score
        ife_em_to_process = set([])
        pdb_id_em_process = set([])
        for ife in ife_no_Q_score:
            pdb_id = ife.split("|")[0]
            if pdb_id in em_pdbs:
                ife_em_to_process.add(ife)
                pdb_id_em_process.add(pdb_id)

        self.logger.info("Found %d pdbs solved by EM but no Q_score" % len(pdb_id_em_process))
        self.logger.info("%s" % sorted(pdb_id_em_process))

        need_to_process = ife_never_processed | ife_em_to_process

        if len(need_to_process) == 0:
            raise core.Skip("No new IFEs to process for CQS")

        return sorted(need_to_process)


    def query(self, session, ife_id):
        if self.fixing:
            # allow re-processing of all IFEs
            return session.query(mod.IfeCqs.ife_id).filter_by(ife_id='nonexistent')
        else:
            # do not re-do IFEs that have already been processed
            return session.query(mod.IfeCqs.ife_id).filter_by(ife_id=ife_id)


    def class_property(self, ifes, name):
        """
        Not sure that this is still used
        """
        return {ife[name] for ife in ifes}


    def ife_info(self, nr_class):
        """
        Not sure that this is still used
        """
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


    def __chain_query__(self, query, info, table=mod.UnitInfo):
        """
        Filtering steps that are used repeatedly in this class.
        info is a dictionary with information about an IFE.
        """
        return query.filter(table.pdb_id == info['pdb']).\
            filter(table.model == info['model']).\
            filter(table.sym_op == info['sym_op']).\
            filter(table.chain.in_(info['chains'])).\
            filter(table.unit_type_id.in_(['rna','dna'])).\
            filter(table.chain_index != None)

    def count_atoms(self, info):
        with self.session() as session:
            query = session.query(mod.UnitCoordinates).\
                join(mod.UnitInfo,
                     mod.UnitInfo.unit_id == mod.UnitCoordinates.unit_id)
            query = self.__chain_query__(query, info)
            counted_atoms = set(['C', 'N', 'O', 'P'])
            count = 0
            for row in query:
                current = 0
                for line in row.coordinates.split('\n'):
                    parts = line.split()
                    if len(parts) >= 2 and parts[2] in counted_atoms:
                        current += 1

                #if not current:
                #    self.logger.error("No atoms in %s" % row.unit_id)
                count += current

            if not count:
                # self.logger.error("No atoms found for %s" % str(info))
                return 100.0

            return float(count)

    def count_clashes(self, info):
        with self.session() as session:
            u1 = aliased(mod.UnitInfo)
            u2 = aliased(mod.UnitInfo)
            query = session.query(mod.UnitClashes).\
                join(u1, u1.unit_id == mod.UnitClashes.unit_id_1).\
                join(u2, u2.unit_id == mod.UnitClashes.unit_id_2).\
                filter(~mod.UnitClashes.atom_name_1.like('%H%')).\
                filter(~mod.UnitClashes.atom_name_2.like('%H%'))

            query = self.__chain_query__(query, info, table=u1)
            query = self.__chain_query__(query, info, table=u2)

            return float(query.count())

    def percent_clash(self, info):
        clashes = self.count_clashes(info)
        if clashes is None:
            return (False, 100)

        atoms = self.count_atoms(info)
        percent_clash = 100 * clashes / atoms
        return (True, percent_clash)


    def average_Q_score(self, info):
        default_avg_Q_score = None  # do not write a value if there is none
        with self.session() as session:
            query = session.query(mod.UnitQuality.Q_score).\
                join(mod.UnitInfo,
                     mod.UnitInfo.unit_id == mod.UnitQuality.unit_id)

            query = self.__chain_query__(query, info)
            if not query.count():
                return (False, default_avg_Q_score)

            values = [r.Q_score for r in query if r.Q_score is not None]
            if not values:
                return (False, default_avg_Q_score)
            return (True, np.mean(values))


    def average_residue_inclusion(self, info):
        default_avg_residue_inclusion = None
        with self.session() as session:
            query = session.query(mod.UnitQuality.residue_inclusion).\
                join(mod.UnitInfo,
                     mod.UnitInfo.unit_id == mod.UnitQuality.unit_id)

            query = self.__chain_query__(query, info)
            if not query.count():
                return (False, default_avg_residue_inclusion)

            values = [r.residue_inclusion for r in query if r.residue_inclusion is not None]
            if not values:
                return (False, default_avg_residue_inclusion)
            return (True, np.mean(values))


    def average_rsr(self, info):
        default_avg_rsr = 40 # on 2017-10-12, maximum observed value was ~31
        with self.session() as session:
            query = session.query(mod.UnitQuality.real_space_r).\
                join(mod.UnitInfo,
                     mod.UnitInfo.unit_id == mod.UnitQuality.unit_id)

            query = self.__chain_query__(query, info)
            if not query.count():
                return (False, default_avg_rsr)

            values = [r.real_space_r for r in query if r.real_space_r is not None]
            if not values:
                return (False, default_avg_rsr)
            return (True, np.mean(values))

    def average_rscc(self, info):
        default_average_rscc = -1 # minimum possible value for rscc
        with self.session() as session:
            query = session.query(mod.UnitQuality.rscc).\
                join(mod.UnitInfo,
                     mod.UnitInfo.unit_id == mod.UnitQuality.unit_id)

            query = self.__chain_query__(query, info)
            if not query.count():
                return (False, default_average_rscc)

            values = [r.rscc for r in query if r.rscc is not None]
            if not values:
                return (False, default_average_rscc)
            return (True, np.mean(values))

    def resolution(self, info):
        with self.session() as session:
            resolution = session.query(mod.PdbInfo.resolution).\
                filter_by(pdb_id=info['pdb']).\
                one().\
                resolution
            if resolution is None:
                return (False, 100)
            return (True, resolution)

    def rfree(self, info):
        with self.session() as session:
            rfree = session.query(mod.PdbQuality).\
                filter_by(pdb_id=info['pdb']).\
                first()

            if rfree is None or rfree.dcc_rfree is None:
                return (False, 1.0)
            return (True, rfree.dcc_rfree)

    # def member_info(self, member):
    #     with self.session() as session:
    #         info = session.query(mod.IfeInfo.pdb_id.label('pdb'),
    #                              mod.IfeInfo.model).\
    #             filter_by(ife_id=member['id']).\
    #             one()
    #         info = row2dict(info)
    #         info.update(member)

    #     with self.session() as session:
    #         query = session.query(mod.ChainInfo.chain_name,
    #                                 mod.IfeChains.is_structured,
    #                                 ).\
    #             join(mod.IfeChains,
    #                     mod.IfeChains.chain_id == mod.ChainInfo.chain_id).\
    #             filter_by(ife_id=member['id'])

    #         if not query.count():
    #             raise core.InvalidState("Could not find chains for %s" %
    #                                     member)

    #         all_chains = [row2dict(c) for c in query]
    #         valid = op.itemgetter('is_structured')
    #         chains = [c['chain_name'] for c in all_chains if valid(c)]
    #         if not chains:
    #             chains = [c['chain_name'] for c in all_chains]

    #     info['chains'] = chains
    #     loader = self._create(IfeLoader)
    #     info['sym_op'] = loader.sym_op(info['pdb'])

    #     return info


    def observed_length_sym_op(self, info):
        """
        Query the unit_info table to determine the correct symmetry
        operator and count the number of observed nucleotides in the first chain

        Until 2025-03-12, this query got distinct chain_index values, which
        meant that chain_index 2 would only get counted once, even if it occurred
        in three chains in the IFE.  Very strange way of counting.
        Now we get distinct chain_index and chain values.
        """
        with self.session() as session:
            query = session.query(mod.UnitInfo.chain_index,mod.UnitInfo.sym_op,mod.UnitInfo.chain).\
                distinct().\
                filter(mod.UnitInfo.pdb_id == info['pdb']).\
                filter(mod.UnitInfo.model == info['model']).\
                filter(mod.UnitInfo.chain.in_(info['chains'])).\
                filter(mod.UnitInfo.unit_type_id.in_(['dna','rna'])).\
                filter(mod.UnitInfo.chain_index != None)

            sym_op_to_count = {}
            for row in query:
                sym_op = row.sym_op
                if not sym_op in sym_op_to_count:
                    sym_op_to_count[sym_op] = 0
                sym_op_to_count[sym_op] += 1

            if '1_555' in sym_op_to_count or len(sym_op_to_count) == 0:
                sym_op = '1_555'
            else:
                sym_op = max(sym_op_to_count.items(), key=op.itemgetter(1))[0]

            nt_count = 0
            for row in query:
                if row.sym_op == sym_op:
                    nt_count += 1

            return nt_count, sym_op

        # until 2024-08-18 instead of unit_type_id the query above used a simpler thing:
        # filter(mod.UnitInfo.unit.in_(['A', 'C', 'G', 'U'])).\
        # That recorded 0 for DNA chains, and ignored modified RNA nucleotides.

    def data(self, ife_id, **kwargs):
        """
        Process one ife id to get data needed to compute the
        composite quality score.

        Parameters
        ----------
        entry : str
            The IFE for which to collect IFE-level composite
            quality score data.

        Returns
        -------
            The required data for the database update step.
        """

        # collect information about this ife in a dictionary
        info = {}
        info['ife_id'] = ife_id

        fields = ife_id.split('|')
        if len(fields) < 3:
            raise core.InvalidState("Invalid ife_id: %s" % ife_id)

        info['pdb'] = fields[0]
        info['model'] = fields[1]

        # use the chains in the ife id, not in accompanying chains
        chains = []
        for chain in ife_id.split("+"):
            fields = chain.split("|")
            if len(fields) >= 3:
                chains.append(fields[2])
        info['chains'] = chains

        # count the number of observed nucleotides in the first chain
        # also get the symmetry operator to use
        observed_length, sym_op = self.observed_length_sym_op(info)
        info['observed_length'] = observed_length
        info['sym_op'] = sym_op

        # calculate components of composite quality score
        # first returned component is Boolean indicating if the data
        # was present or not, but that is not actually used
        a, info['average_rsr'] = self.average_rsr(info)
        b, info['average_rscc'] = self.average_rscc(info)
        c, info['percent_clash'] = self.percent_clash(info)
        d, info['rfree'] = self.rfree(info)
        e, info['resolution'] = self.resolution(info)
        f, info['average_Q_score'] = self.average_Q_score(info)
        g, info['average_residue_inclusion'] = self.average_residue_inclusion(info)

        if False and (self.testing or self.fixing):
            # print('Processed %s' % ife_id)
            # load the data from the database if it is there already
            with self.session() as session:
                query = session.query(mod.IfeCqs).filter(mod.IfeCqs.ife_id==ife_id)

                if query.count() > 1:
                    self.logger.info('Error: multiple entries for %s' % ife_id)
                if query.count() == 1:
                    for row in query:
                        if abs(row.average_rsr - info['average_rsr']) > 0.01:
                            self.logger.info('Change in %s: old rsr: %s, new rsr: %s' % (row.ife_id,row.average_rsr, info['average_rsr']))
                        if abs(row.average_rscc - info['average_rscc']) > 0.01:
                            self.logger.info('Change in %s: old rscc: %s, new rscc: %s' % (row.ife_id,row.average_rscc, info['average_rscc']))
                        if abs(row.percent_clash - info['percent_clash']) > 0.01:
                            self.logger.info('Change in %s: old clash: %s, new clash: %s' % (row.ife_id,row.percent_clash, info['percent_clash']))
                        if abs(row.rfree - info['rfree']) > 0.01:
                            self.logger.info('Change in %s: old rfree: %s, new rfree: %s' % (row.ife_id,row.rfree, info['rfree']))
                        if abs(row.resolution - info['resolution']) > 0.01:
                            self.logger.info('Change in %s: old resolution: %s, new resolution: %s' % (row.ife_id,row.resolution, info['resolution']))
                        if row.obs_length < info['observed_length']:
                            self.logger.info('Change in %s: longer: old observed_length: %s, new observed_length: %s' % (row.ife_id,row.obs_length, info['observed_length']))
                        if row.obs_length > info['observed_length']:
                            self.logger.info('Change in %s: shorter: old observed_length: %s, new observed_length: %s' % (row.ife_id,row.obs_length, info['observed_length']))

        if info['average_Q_score'] and info['average_residue_inclusion']:
            self.logger.info("%s has average_Q_score %10.3f and average_residue_inclusion %10.3f" % (info['ife_id'],info['average_Q_score'],info['average_residue_inclusion']))
        else:
            self.logger.info("%s has average_Q_score %s and average_residue_inclusion %s" % (info['ife_id'],info['average_Q_score'],info['average_residue_inclusion']))

        if self.testing:
            print("Processed %s for testing" % ife_id)
            print("%s has average_Q_score %s and average_residue_inclusion %s" % (info['ife_id'],info['average_Q_score'],info['average_residue_inclusion']))
            self.logger.info("info values %s" % info)
            raise core.Skip('Not writing ife_cqs entries yet')
        elif self.fixing:
            print("%s has average_Q_score %s and average_residue_inclusion %s" % (info['ife_id'],info['average_Q_score'],info['average_residue_inclusion']))
            # modify the line below for whatever you are fixing
            yield mod.IfeCqs(
                ife_id = info['ife_id'],
                average_Q_score = info['average_Q_score'],
                average_residue_inclusion = info['average_residue_inclusion'])
        else:
            # production situation
            yield mod.IfeCqs(
                ife_id = info['ife_id'],
                obs_length = info['observed_length'],
                clashscore = 100,                        # default value for some reason
                average_rsr = info['average_rsr'],
                average_rscc = info['average_rscc'],
                percent_clash = info['percent_clash'],
                rfree = info['rfree'],
                resolution = info['resolution'],
                average_Q_score = info['average_Q_score'],
                average_residue_inclusion = info['average_residue_inclusion'])
