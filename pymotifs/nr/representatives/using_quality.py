import abc

import numpy as np

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import row2dict

from pymotifs.constants import MANUAL_IFE_REPRESENTATIVES
from pymotifs.constants import WORSE_THAN_MANUAL_IFE_REPRESENTATIVES

from pymotifs.representatives.core import Representative


class QualityBase(Representative):
    """
    Find representatives using quality metrics provided by PDB.
    """
    hardcoded = MANUAL_IFE_REPRESENTATIVES
    worse = WORSE_THAN_MANUAL_IFE_REPRESENTATIVES

    __metaclass__ = abc.ABCMeta()

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
        current = members[0]
        hardcoded = self.find_hardcoded(members)

        if not hardcoded:
            self.logger.debug("No hardcoded representative to use")
            return list(members)

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
        self.logger.info("Selecting representative for %s", given)
        with_quailty = self.load_quality(given['members'])
        candidates = self.select_candidates(with_quailty)
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
    name = 'compscore'

    def has_quality(self, member):
        return 'real_space_r' in member['has'] and \
            'rscc' in member['has'] and \
            'pdb_rfree' in member['has'] and \
            'resolution' in member['has']

    def as_quality(self, entries):
        def avg_of(name):
            return sum(e[name] for e in entries if e[name]) / len(entries)

        def has_entry(name):
            return any(e[name] is not None for e in entries)

        def percent_of(name):
            return 100 * sum(e[name] for e in entries) / len(entries)

        def first_value(name):
            return entries[0][name]

        def assign(name, default, function, tracking):
            if has_entry(name):
                tracking.add(name)
                return (function(name), tracking)
            return (default, tracking)

        resolution, has = assign('resolution', 100, first_value, set())
        rfree, has = assign('rfree', 1, first_value, has)
        average_rsr, has = assign('real_space_r', 1, avg_of, has)
        average_rscc, has = assign('rscc', 0, avg_of, has)

        percent_clash = 100
        if entries[0]['clash_score'] is not None:
            percent_clash, has = assign('clash_count', 0, percent_of, has)

        return {
            'resolution': resolution,
            'percent_clash': percent_clash,
            'average_rsr': average_rsr,
            'average_rscc': average_rscc,
            'rfree': rfree,
            'has': has,
        }

    def load_quality(self, members):
        for member in members:
            chains = [mem['chain'] for mem in members if members['structured']]
            if not chains:
                chains = [mem['chain'] for mem in members]

            with self.session() as session:
                query = session.query(
                    mod.UnitQuality.real_space_r,
                    mod.UnitQuality.clash_count,
                    mod.UnitQuality.rscc,
                    mod.PdbQuality.dcc_rfree.label('rfree'),
                    mod.PdbQuality.clashscore,
                    mod.PdbInfo.resolution,
                ).join(mod.UnitInfo
                       mod.UnitQuality.unit_id == mod.UnitInfo.unit_id).\
                    join(mod.PdbInfo, mod.PdbInfo.pdb_id == mod.UnitInfo.pdb_id).\
                    filter(mod.UnitInfo.pdb_id == member['pdb']).\
                    filter(mod.UnitInfo.model == member['model']).\
                    filter(mod.UnitInfo.sym_op == member['sym_op']).\
                    filter(mod.UnitInfo.chain.in_(chains)).\
                    filter(mod.UnitInfo.unit.in_(['A', 'C', 'G', 'U']))

            member['quality'] = self.as_quality([row2dict(r) for r in query])
        return member

    def compscore(self, member):
        average = np.mean(member['resolution'],
                          member['percent_clash'],
                          10 * member['average_rsr'],
                          10 * (1 - member['average_rscc']),
                          10 * member['rfree'])

        return 100 * average

    def sort_by_quality(self, members):
        return sorted(members, key=self.compscore, reverse=True)
