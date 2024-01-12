import abc
import operator as op

import numpy as np

from sqlalchemy import func
from sqlalchemy.orm import aliased

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import row2dict

from pymotifs.constants import COMPSCORE_COEFFICENTS
from pymotifs.constants import MANUAL_IFE_REPRESENTATIVES
from pymotifs.constants import WORSE_THAN_MANUAL_IFE_REPRESENTATIVES

from pymotifs.ife.helpers import IfeLoader

from .core import Representative

import inspect
import traceback


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
        already_ordered = {m['id'] for m in ordered}
        rest = [m for m in given if m['id'] not in already_ordered]
        return ordered + self.sort_by_quality(rest)

    def __call__(self, given):
        # # Get information about the __call__ function
        # caller_frame = inspect.currentframe().f_back
        # caller_function_name = caller_frame.f_code.co_name
        # caller_module = inspect.getmodule(caller_frame).__name__

        # # Use traceback to get the calling code line
        # stack = traceback.extract_stack()
        # calling_code_line = stack[-2][1]
        # self.logger.info("__call__ method called from function '%s' in module '%s' at line %s" % (caller_function_name, caller_module, calling_code_line))
        # self.logger.info("show the runnng message 1: %s", stack[-3])
        # self.logger.info("show the runnng message 2: %s", stack[-2])
        # self.logger.info("show the runnng message 3: %s", stack[-1])
        # self.logger.info("show the runnng message 4: %s", stack[0])
        # self.logger.info("show the runnng message 5: %s", stack[1])
        # self.logger.info("show the runnng message 6: %s", stack[2])


        # self.logger.info("Selecting representative for %s", given['name']['full'])
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
    """
    method = 'compscore'

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
                self.logger.error("No atoms found for %s" % str(info))
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

    def has_quality(self, member):
        has_quality = member['quality']['has']
        return 'real_space_r' in has_quality and \
            'rscc' in has_quality and \
            'rfree' in has_quality and \
            'resolution' in has_quality

    def select_candidates(self, members):
        return [m for m in members if self.has_quality(m)]

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
        experimental = float(info['max_length'])
        #return (True, (observed - 1) / experimental)
        self.logger.debug("info: %s" % info)
        self.logger.debug("observed: %s" % observed)
        self.logger.debug("experimental: %s" % experimental)
        return (True, (1 - (observed / experimental)))

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

    def experimental_length(self, members):
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

    def load_quality(self, members):
        """
        This will load and store all quality data for the given list of members
        of the EC.
        """

        parameters = [
            'resolution',
            'percent_clash',
            'average_rsr',
            'average_rscc',
            'rfree',
            'fraction_unobserved',
        ]

        experimental_length = self.experimental_length(members)

        for member in members:
            info = self.member_info(member)
            info['max_length'] = experimental_length
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

            #if has:
                #data['has'].add(('max_length', experimental_length))

            data['max_length'] = experimental_length
            data['obs_length'] = info['length']

            member['quality'] = data
        return members

    def compscore(self, member):
        """
        Compute composite quality score using six indicators weighted by various coefficients
        set in constants.py
        In the future, we may just use the quality instead of doing the database query ___ 9/1/2023
        """
        quality = member['quality']
        ife_id = member['id']
        nr_class = member['name']

        with self.session() as session:
            query = session.query(mod.NrCqs.ife_id,
                                  mod.NrCqs.nr_name,
                                  mod.NrCqs.composite_quality_score.label('cqs'),
                                  ).\
                filter(mod.NrCqs.ife_id == ife_id).\
                filter(mod.NrCqs.nr_name == nr_class).\
                limit(1)

            for result in query:
                cqs = result.cqs
        #self.logger.info("compscore function of using_quality")
        try:
            cqs
        except NameError:
            pass
            cqs = COMPSCORE_COEFFICENTS['resolution'] * quality['resolution']
            cqs += COMPSCORE_COEFFICENTS['percent_clash'] * quality['percent_clash']
            cqs += COMPSCORE_COEFFICENTS['average_rsr'] * quality['average_rsr']
            cqs += COMPSCORE_COEFFICENTS['average_rscc'] * (1 - quality['average_rscc'])
            cqs += COMPSCORE_COEFFICENTS['rfree'] * quality['rfree']
            cqs += COMPSCORE_COEFFICENTS['fraction_unobserved'] * quality['fraction_unobserved']

        self.logger.debug("test compscore: %s" % cqs)

        if cqs < 0:
            raise core.InvalidState("Invalid compscore (%s) for %s" % (cqs, member))

        return cqs


    def sort_by_quality(self, members):
        return sorted(members, key=self.compscore, reverse=False)
        #return sorted(members, key=self.compscore, reverse=True)

    def __chain_query__(self, query, info, table=mod.UnitInfo):
        return query.filter(table.pdb_id == info['pdb']).\
            filter(table.model == info['model']).\
            filter(table.sym_op == info['sym_op']).\
            filter(table.chain.in_(info['chains'])).\
            filter(table.unit.in_(['A', 'C', 'G', 'U'])).\
            filter(table.chain_index != None)
