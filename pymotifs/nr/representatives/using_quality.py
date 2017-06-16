import abc

import numpy as np

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import row2dict

from pymotifs.constants import MANUAL_IFE_REPRESENTATIVES
from pymotifs.constants import WORSE_THAN_MANUAL_IFE_REPRESENTATIVES

from pymotifs.ife.helpers import IfeLoader

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
    method = 'compscore'

    def has_quality(self, member):
        has_quality = member['quality']['has']
        return 'real_space_r' in has_quality and \
            'rscc' in has_quality and \
            'rfree' in has_quality and \
            'resolution' in has_quality

    def select_candidates(self, members):
        return [m for m in members if self.has_quality(m)]

    def as_quality(self, atoms, entries):
        if not entries:
            raise core.InvalidState("No entries to compute quality for")

        def avg_of(name):
            return np.mean([e[name] for e in entries])

        def has_entry(name):
            return any(e[name] is not None for e in entries)

        def first_value(name):
            values = set(entry[name] for entry in entries)
            if len(values) > 1:
                raise core.InvalidState("Excpected only a single value for %s"
                                        % name)
            return values.pop()

        def assign(name, default, function, tracking):
            if has_entry(name):
                tracking.add(name)
                return (function(name), tracking)
            return (default, tracking)

        resolution, has = assign('resolution', 100, first_value, set())
        rfree, has = assign('rfree', 1, first_value, has)
        average_rsr, has = assign('real_space_r', 1, avg_of, has)
        average_rscc, has = assign('rscc', 0, avg_of, has)

        clash_count = 0.0
        percent_clash = 100
        if entries[0]['clashscore'] is not None:
            clash_count = sum(e['clash_count'] for e in entries)
            percent_clash = clash_count / atoms
            has.add('clashscore')

        return {
            'resolution': resolution,
            'clash': clash_count,
            'percent_clash': percent_clash,
            'average_rsr': average_rsr,
            'average_rscc': average_rscc,
            'rfree': rfree,
            'has': has,
            'atoms': atoms,
        }

    def member_info(self, member):
        with self.session() as session:
            info = session.query(mod.IfeInfo.pdb_id.label('pdb'),
                                 mod.IfeInfo.model).\
                filter_by(ife_id=member['id']).\
                one()
            info = row2dict(info)

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
                chains = [c['chain_name'] for c in all_chains if c['is_structured']]
                if not chains:
                    chains = [c['chain_name'] for c in all_chains]

            info['chains'] = chains

            loader = self._create(IfeLoader)
            info['sym_op'] = loader.sym_op(info['pdb'])

            return info

    def __chain_query__(self, query, info):
        return query.filter(mod.UnitInfo.pdb_id == info['pdb']).\
            filter(mod.UnitInfo.model == info['model']).\
            filter(mod.UnitInfo.sym_op == info['sym_op']).\
            filter(mod.UnitInfo.chain.in_(info['chains'])).\
            filter(mod.UnitInfo.unit.in_(['A', 'C', 'G', 'U']))

    def count_atoms(self, info):
        with self.session() as session:
            query = session.query(mod.UnitCoordinates).\
                join(mod.UnitInfo,
                     mod.UnitInfo.unit_id == mod.UnitCoordinates.unit_id)
            query = self.__chain_query__(query, info)
            counted_atoms = set(['C', 'N', 'O'])
            count = 0
            for row in query:
                current = 0
                for line in row.coordinates.split('\n'):
                    parts = line.split()
                    if parts[2] in counted_atoms:
                        current += 1

                if not current:
                    raise core.InvalidState("No atoms in %s" % row.unit_id)
                count += current

            if not count:
                raise core.InvalidState("No atoms found for %s" % str(info))

            return float(count)

    def load_quality(self, members):
        for member in members:
            info = self.member_info(member)
            atoms = self.count_atoms(info)

            with self.session() as session:
                query = session.query(
                    mod.UnitQuality.real_space_r,
                    mod.UnitQuality.clash_count,
                    mod.UnitQuality.rscc,
                    mod.PdbQuality.dcc_rfree.label('rfree'),
                    mod.PdbQuality.clashscore,
                    mod.PdbInfo.resolution,
                ).join(mod.UnitInfo,
                       mod.UnitQuality.unit_id == mod.UnitInfo.unit_id).\
                    join(mod.PdbInfo,
                         mod.PdbInfo.pdb_id == mod.UnitInfo.pdb_id).\
                    join(mod.PdbQuality,
                         mod.PdbQuality.pdb_id == mod.PdbInfo.pdb_id)
                query = self.__chain_query__(query, info)

                entries = [row2dict(r) for r in query]
            member['quality'] = self.as_quality(atoms, entries)
        return members

    def compscore(self, member):
        quality = member['quality']
        average = np.mean([quality['resolution'],
                           quality['percent_clash'],
                           10 * quality['average_rsr'],
                           10 * (1 - quality['average_rscc']),
                           10 * quality['rfree']])

        return 100 * average

    def sort_by_quality(self, members):
        return sorted(members, key=self.compscore, reverse=True)
