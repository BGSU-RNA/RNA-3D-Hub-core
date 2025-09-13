import abc
# import numpy as np
# import operator as op

from sqlalchemy import func
# from sqlalchemy.orm import aliased

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import row2dict

from pymotifs.constants import COMPSCORE_COEFFICENTS
from pymotifs.constants import CQS2_COEFFICIENTS
from pymotifs.constants import MANUAL_IFE_REPRESENTATIVES
from pymotifs.constants import WORSE_THAN_MANUAL_IFE_REPRESENTATIVES

# from pymotifs.ife.helpers import IfeLoader

from .core import Representative


class QualityBase(Representative):
    """
    Find representatives using quality metrics provided by PDB.
    """
    hardcoded = MANUAL_IFE_REPRESENTATIVES
    worse = WORSE_THAN_MANUAL_IFE_REPRESENTATIVES

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def has_quality(self, member):
        pass

    @abc.abstractmethod
    def load_quality(self, members):
        """
        Load all required quality data for the given members.
        """
        pass

    @abc.abstractmethod
    def sort_by_quality(self, members):
        """
        This must provide a linear ordering to all members in the group. The
        one which is listed first will be used as the representative. In
        addition it should be tolerant of missing data and the like as it will
        be called twice. First with X-ray structures with all quality data and
        then later with non-xray structures that may be missing data.
        """
        pass

    @abc.abstractmethod
    def select_candidates(self, members):
        return members

    def find_hardcoded(self, members):
        ids = {m['id'] for m in members}

        found = ids.intersection(self.hardcoded)
        if len(found) > 1:
            self.logger.error("More than one hardcoded, using quality")
            return None

        if not found:
            self.logger.debug("No hardcoded representative")
            return None

        selected = next(m for m in members if m['id']in found)
        self.logger.info("Found hardcoded representative %s", selected)
        return selected

    def filter_by_method(self, members):
        def is_good_xray(member):
            return member['method'] == 'X-RAY DIFFRACTION' and \
                member['resolution'] <= 4.0 and \
                self.has_quality(member)

        if any(is_good_xray(m) for m in members):
            return [m for m in members if is_good_xray(m)]
        return list(members)

    def filter_by_nts(self, members):
        best = max(m['length'] for m in members)
        return [m for m in members if m['length'] >= 0.75 * best]

    def filter_by_resolution(self, members):
        best = min(m['resolution'] for m in members)
        return [m for m in members if abs(m['resolution'] - best) <= 0.2]

    def use_hardcoded(self, members):
        hardcoded = self.find_hardcoded(members)
        if not hardcoded:
            self.logger.debug("No hardcoded representative to use")
            return list(members)

        if not members:
            return []

        current = members[0]
        if current == hardcoded:
            self.logger.info("Hardcoded and chosen agree")
            return list(members)

        if current['id'] not in self.worse:
            self.logger.warning("Automatically selected %s not in"
                                " WORSE_THAN_MANUAL_IFE_REPRESENTATIVES",
                                current)

        return self.insert_as_representative(hardcoded, members)

    def final_ordering(self, ordered, given):
        """
        This method is still used.
        """
        already_ordered = {m['id'] for m in ordered}

        rest = [m for m in given if m['id'] not in already_ordered]

        # already_ordered appears to be empty in all cases
        # self.logger.info('final_ordering: already_ordered: %s' % already_ordered)
        # self.logger.info('final_ordering: rest: %s' % rest)

        return ordered + self.sort_by_quality(rest)

    def __call__(self, given):
        """
        This is the abstract method when ordering members of a class.
        For builder.py, this is called when rep_finder is called.
        See below the class with method = 'compscore' for equivalence classes.
        """

        self.logger.info("Ordering members of %s", given['name']['full'])
        with_quality = self.load_quality(given['members'])
        candidates = self.select_candidates(with_quality)
        ordered_by_quality = self.sort_by_quality(candidates)
        with_representative = self.use_hardcoded(ordered_by_quality)
        return self.final_ordering(with_representative, given['members'])


class QualityMetrics(QualityBase):
    method = 'quality-metrics'

    def has_quality(self, member):
        required = set(['rsrz', 'backbone', 'clashscore'])
        return required.issubset(member['quality']['has'])

    def load_quality(self, members):
        def as_quality(data):
            return {
                'has': {key for key, value in data.items() if value},
                'rsrz': data.get('rsrz') or 100,
                'backbone': data.get('backbone') or 100,
                'clashscore': data.get('clashscore') or 500,
            }

        known = {m['pdb'] for m in members}
        with self.session() as session:
            query = session.query(mod.PdbQuality.pdb_id,
                                  mod.PdbQuality.percent_rsrz_outliers.
                                  label('rsrz'),
                                  mod.PdbQuality.clashscore,
                                  mod.PdbQuality.percent_rota_outliers.
                                  label('backbone'),
                                  ).\
                filter(mod.PdbQuality.pdb_id.in_(known))

            measures = {}
            for result in query:
                result = row2dict(result)
                pdb_id = result.pop('pdb_id')
                measures[pdb_id] = as_quality(result)

        for member in members:
            pdb_id = member['pdb']
            member['quality'] = measures.get(pdb_id, as_quality({}))

        return members

    def sort_by_quality(self, members):
        def key(member):
            quality = member['quality']
            quality_factor = quality['rsrz'] * \
                quality['clashscore'] * \
                quality['backbone'] * \
                pow(member['resolution'], 4)
            size_factor = member['length'] + member['bp']
            return quality_factor / size_factor

        return sorted(members, key=key)

    def select_candidates(self, members):
        best_method = self.filter_by_method(members)
        best_nts = self.filter_by_nts(best_method)
        best_resolution = self.filter_by_resolution(best_nts)
        if not best_resolution:
            raise core.InvalidState("Nothing with good resolution")
        return best_resolution


class CompScore(QualityBase):
    """
    This implements a composite scoring metric that should provide a good
    linear ordering of IFE's.

    Apparently several of the methods here are no longer used, having
    been replaced by the nr.cqs stage that stores the calculations in the
    nr_cqs table, which is faster.

    Or ... maybe replaced by the ife/cqs stage which stores in ife_cqs table.

    Unfortunately, that makes it harder to follow this code.
    """
    method = 'compscore'

    """
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
            if query.count() < 0:
                raise core.InvalidState("Negative clashes: %s" % info)
            return float(query.count())

    def percent_clash(self, info):
        clashes = self.count_clashes(info)
        if clashes is None:
            return (False, 100)

        atoms = self.count_atoms(info)
        percent_clash = 100 * clashes / atoms
        return (True, percent_clash)

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
    """

    def observed_length(self, info):
        self.logger.debug("info: %s" % info)
        with self.session() as session:
            #query = session.query(mod.UnitInfo.unit_id).\
            query = session.query(mod.UnitInfo.chain_index).\
                distinct()
            query = self.__chain_query__(query, info)
            return query.count()

    def fraction_unobserved(self, info):
        observed = float(self.observed_length(info))
        experimental = float(info['maximum_experimental_length'])

        return (True, (1 - (observed / experimental)))

    def has_quality(self, member):
        """
        """
        has_quality = member['quality']['has']
        return 'real_space_r' in has_quality and \
            'rscc' in has_quality and \
            'rfree' in has_quality and \
            'resolution' in has_quality

    def select_candidates(self, members):
        """
        This method is still used.
        """
        return [m for m in members if "composite_quality_score" in m]

    def select_candidates_old(self, members):
        """
        """
        return [m for m in members if self.has_quality(m)]

    def member_info(self, member):
        """
        """

        ife_id = member['id']
        fields = ife_id.split('|')
        if len(fields) < 3:
            raise core.InvalidState("Invalid ife_id: %s" % ife_id)

        pdb_id = fields[0]
        model = fields[1]
        member['pdb'] = pdb_id
        member['model'] = model

        # no need to query to get this information
        # with self.session() as session:
        #     info = session.query(mod.IfeInfo.pdb_id.label('pdb'),
        #                          mod.IfeInfo.model).\
        #                     filter_by(ife_id=member['id']).\
        #                     filter(mod.IfeInfo.new_style == True).\
        #                     one()
        #     info = row2dict(info)
        #     info.update(member)

        # find out which chains in this ife are structured
        # but ... that cannot be important for IFE ranking, because many
        # chains are not structured
        # Just trust the chains that are in the ife_id
        # with self.session() as session:
        #     query = session.query(mod.ChainInfo.chain_name,
        #                             mod.IfeChains.is_structured,
        #                             ).\
        #         join(mod.IfeChains,
        #                 mod.IfeChains.chain_id == mod.ChainInfo.chain_id).\
        #         filter_by(ife_id=member['id'])

        #     if not query.count():
        #         raise core.InvalidState("Could not find chains for %s" %
        #                                 member)

        #     all_chains = [row2dict(c) for c in query]
        #     valid = op.itemgetter('is_structured')
        #     chains = [c['chain_name'] for c in all_chains if valid(c)]
        #     if not chains:
        #         chains = [c['chain_name'] for c in all_chains]

        chains = []
        for chain in ife_id.split("+"):
            fields = chain.split("|")
            if len(fields) >= 3:
                chains.append(fields[2])

        # loader = self._create(IfeLoader)
        # info['sym_op'] = loader.sym_op(info['pdb'])

        return member

    def maximum_experimental_length(self, members):
        """
        Turns out that maximum experimental length is not stable enough
        because some chains have the incorrect experimental sequence
        reported for them, much longer than the actual molecule used in
        the experiment.  Better to use maximum number of observed nucleotides.
        """
        ids = [m['id'] for m in members]
        with self.session() as session:
            query = session.query(func.sum(mod.ExpSeqInfo.length).label('length')).\
                join(mod.ExpSeqChainMapping,
                     mod.ExpSeqChainMapping.exp_seq_id == mod.ExpSeqInfo.exp_seq_id).\
                join(mod.IfeChains,
                     mod.IfeChains.chain_id == mod.ExpSeqChainMapping.chain_id).\
                filter(mod.IfeChains.ife_id.in_(ids)).\
                group_by(mod.IfeChains.ife_id)
            return max(r.length for r in query)

    def load_quality_old(self, members):
        """
        Load and store all quality data for the given list of members of the EC.
        """

        parameters = [
            'resolution',
            'percent_clash',
            'average_rsr',
            'average_rscc',
            'rfree',
            'fraction_unobserved',
        ]

        # get the maximum experimental length across this equivalence class
        maximum_experimental_length = self.maximum_experimental_length(members)

        for member in members:
            # query database for pdb id, model, symmetry operator
            info = self.member_info(member)

            info['maximum_experimental_length'] = maximum_experimental_length
            data = {'has': set()}
            for index, name in enumerate(parameters):
                method = getattr(self, name)
                has, value = method(info)
                if value < 0 and has and name != 'average_rscc':
                    raise core.InvalidState("%s should be positive: %s %s" %
                                            (name, value, info))

                data[name] = value
                if has:
                    data['has'].add((name, index))

            data['maximum_experimental_length'] = maximum_experimental_length
            data['obs_length'] = info['length']

            member['quality'] = data

        return members

    def load_quality(self, members):
        """
        Developed as part of nr/cqs.py, but moved here where it belongs.
        That code has been used for a long time to populate the nr_cqs table.
        cqs needs to be calculated here, when an equivalence class is first created.
        """

        ife_list = [m['id'] for m in members]

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
                mod.IfeCqs.average_Q_score,
                mod.IfeCqs.average_residue_inclusion,
                mod.PdbInfo.experimental_technique
                ).\
                join(mod.IfeInfo, mod.IfeInfo.ife_id == mod.IfeCqs.ife_id).\
                join(mod.PdbInfo, mod.PdbInfo.pdb_id == mod.IfeInfo.pdb_id).\
                filter(mod.IfeCqs.ife_id.in_(ife_list))

            ife_id_to_data = {}
            for result in query:
                entry = row2dict(result)
                ife_id = entry['ife_id']
                ife_id_to_data[ife_id] = entry

        # find the maximum number of observed nucleotides in this equivalence class
        max_observed = 0
        for ife_id in ife_list:
            if ife_id in ife_id_to_data:
                obs_length = ife_id_to_data[ife_id]['obs_length']
                if obs_length > max_observed:
                    max_observed = obs_length

        # calculate cqs for each member of the equivalence class
        for member in members:
            ife_id = member['id']
            if ife_id in ife_id_to_data:
                obs_length = ife_id_to_data[ife_id]['obs_length']

                if max_observed == 0:
                    fraction_unobserved = 1
                else:
                    fraction_unobserved = 1 - (obs_length / max_observed)

                member['max_observed'] = max_observed
                member['fraction_unobserved'] = fraction_unobserved
                member['percent_observed'] = 1 - fraction_unobserved

                # calculate original composite_quality_score
                compscore =  COMPSCORE_COEFFICENTS['resolution'] * ife_id_to_data[ife_id]['resolution']
                compscore += COMPSCORE_COEFFICENTS['percent_clash'] * ife_id_to_data[ife_id]['percent_clash']
                compscore += COMPSCORE_COEFFICENTS['fraction_unobserved'] * fraction_unobserved
                compscore += COMPSCORE_COEFFICENTS['average_rsr'] * ife_id_to_data[ife_id]['average_rsr']
                compscore += COMPSCORE_COEFFICENTS['average_rscc'] * (1 - ife_id_to_data[ife_id]['average_rscc'])
                compscore += COMPSCORE_COEFFICENTS['rfree'] * ife_id_to_data[ife_id]['rfree']

                member['composite_quality_score'] = compscore

                # calculate cqs2 using the same code as in ife_fixer.py
                ife_data = ife_id_to_data[ife_id]
                cqs2 = CQS2_COEFFICIENTS['resolution'] * ife_data['resolution']
                cqs2 += CQS2_COEFFICIENTS['percent_clash'] * ife_data['percent_clash']
                cqs2 += CQS2_COEFFICIENTS['fraction_unobserved'] * fraction_unobserved

                if ife_data['experimental_technique'] in ['ELECTRON MICROSCOPY']:
                    if ife_data['average_Q_score'] is None:
                        # match half of the penalty given to non-EM structures missing key data
                        cqs2 += 150
                    else:
                        cqs2 += CQS2_COEFFICIENTS['average_Q_score'] * (1-ife_data['average_Q_score'])
                    if ife_data['average_residue_inclusion'] is None:
                        # match half of the penalty given to non-EM structures missing key data
                        cqs2 += 141
                    else:
                        cqs2 += CQS2_COEFFICIENTS['average_residue_inclusion'] * (1-ife_data['average_residue_inclusion'])
                    cqs2 += CQS2_COEFFICIENTS['constant']
                else:
                    # when these are missing from an x-ray or NMR structure, the addition to the
                    # cqs2 is 6.4*40 + 7*2 + 21*1 = 291
                    cqs2 += CQS2_COEFFICIENTS['average_rsr'] * ife_data['average_rsr']
                    cqs2 += CQS2_COEFFICIENTS['average_rscc'] * (1 - ife_data['average_rscc'])
                    cqs2 += CQS2_COEFFICIENTS['rfree'] * ife_data['rfree']

                member['cqs2'] = cqs2

            else:
                self.logger.info("No data in ife_cqs query for %s" % ife_id)
                member['composite_quality_score'] = 999

        return members


    def sort_by_quality(self, members):
        """
        Sort by computed composite quality score
        """
        return sorted(members, key=lambda x: x['cqs2'])
        return sorted(members, key=lambda x: x['composite_quality_score'])

    def sort_by_quality_old(self, members):
        """
        This is the method that calls compscore.
        """
        return sorted(members, key=self.compscore, reverse=False)

    def compscore(self, member):
        """
        Does not seem that this code gets run.

        Look up the quality score from the nr_cqs table.
        The value of CQS is computed based on fraction observed, which is
        computed based on the largest sequence length in the class.
        So it needs the nr_name available in order to retrieve the CQS.

        Not sure that this function ever gets called, because the
        structure of the code now is that the data for nr_cqs is computed
        at the same time that the ranking is being computed for the first time.
        """

        quality = member['quality']
        ife_id = member['id']
        nr_name = member['nr_name']

        cqs = 999  # default value larger than would ever happen

        with self.session() as session:
            query = session.query(mod.NrCqs.composite_quality_score.label('cqs')).\
                filter(mod.NrCqs.ife_id == ife_id).\
                filter(mod.NrCqs.nr_name == nr_name).\
                limit(1)

            for result in query:
                cqs = result.cqs

        if cqs == 999:
            # sometimes nr_class is equal to the last chain name, then of course this won't work
            self.logger.info('Not able to retrieve CQS for %s in %s' % (ife_id, nr_name))
            pass
            cqs = COMPSCORE_COEFFICENTS['resolution'] * quality['resolution']
            cqs += COMPSCORE_COEFFICENTS['percent_clash'] * quality['percent_clash']
            cqs += COMPSCORE_COEFFICENTS['average_rsr'] * quality['average_rsr']
            cqs += COMPSCORE_COEFFICENTS['average_rscc'] * (1 - quality['average_rscc'])
            cqs += COMPSCORE_COEFFICENTS['rfree'] * quality['rfree']
            cqs += COMPSCORE_COEFFICENTS['fraction_unobserved'] * quality['fraction_unobserved']
        elif cqs < 0:
            raise core.InvalidState("Invalid compscore (%s) for %s" % (cqs, member))
        else:
            self.logger.info('Successfully retrieved CQS %0.2f for %s in %s' % (cqs, ife_id, nr_name))

        return cqs

