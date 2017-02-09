import abc
import copy
import inspect
import operator as op
import itertools as it

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import row2dict
from pymotifs.utils import known_subclasses

from pymotifs.constants import NR_BP_PERCENT_INCREASE
from pymotifs.constants import NR_LENGTH_PERCENT_INCREASE
from pymotifs.constants import NR_ALLOWED_METHODS
from pymotifs.constants import MANUAL_IFE_REPRESENTATIVES
from pymotifs.constants import WORSE_THAN_MANUAL_IFE_REPRESENTATIVES


def bp_per_nt(chain):
    """Function to use for sorting by bps/nt. Deals with things where there
    are 0 bps or nts. It must have bps, length and id entry.

    :param dict chain: Chain to compute a sorting key for.
    :returns: A tuple that can be used to sort the chains.
    """
    ratio = 0
    if chain['bp'] and chain['length']:
        ratio = float(chain['bp']) / float(chain['length'])
    resolution = chain.get('resolution')
    if resolution:
        resolution = resolution * -1
    return (ratio, resolution, chain['id'])


class Representative(core.Base):
    """Just a base class for all things that find representatives. This is only
    for book keeping purposes and implements nothing.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def method(self):
        """The method name for this representative finder. This should be a
        unique string among all representative subclasses.
        """
        return None

    def filter_group_by_method(self, group, methods=NR_ALLOWED_METHODS):
        """This will filter the group that is being examiend to just a copy of
        the parent, and members entries. This is done so that we have a group
        that we can manipulate without messing up parts elsewhere. In addition,
        the members of the group will be filtered to only those with allowed
        methods. These are the methods listed in the methods set.

        Parameters
        ----------
        group : dict
            A group dictonary that must have 'parent' and 'members' entries.
        methods : set
            A set of method names that are allowed. If no members have the
            given method then all are used.

        Returns
        -------
        copied : dict
            A group dictonary with only the 'parent' and 'members' entries.
        """

        meth = op.itemgetter('method')
        members = [ife for ife in group['members'] if meth(ife) in methods]
        if not members:
            members = group['members']

        return {
            'parent': copy.deepcopy(group.get('parent', [])),
            'members': copy.deepcopy(members),
        }

    def insert_as_representative(self, representative, members, sort=None):
        ordered = [representative]
        to_add = [m for m in members if m['id'] != representative['id']]
        if sort:
            to_add.sort(key=sort, reverse=True)
        ordered.extend(to_add)
        return ordered


def known():
    finders = []
    for subclass in known_subclasses(Representative, globals()):
        finders.append((subclass.method, subclass))
    return finders


def fetch(name):
    """Get the class that implements the given method name.

    Parameters
    ----------
    name : str
        Method name to use

    Returns
    -------
    finder : class
        A class that implements the given method.
    """
    for key, value in globals().items():
        if inspect.isclass(value) and issubclass(value, Representative) and \
                getattr(value, 'method', None) == name:
            return value
    raise ValueError("Unknown method %s" % name)


class Naive(Representative):
    """This will select the naive representative. This just selects the chain
    with the most bp/nt and uses that as the representative.
    """
    method = 'naive'

    def __call__(self, group):
        return sorted(group['members'], key=bp_per_nt, reverse=True)


class Increase(Representative):
    """A class to find the representative for a group of ifes. This will find
    the best in terms of bp/nt and attempt to find all those with more bp's and
    nts in the set.
    """
    method = 'percent-increase'

    def initial_representative(self, group):
        """Find the best chain terms of bps/nts. This is the starting point for
        finding the representative in a set of ifes. This method is naive
        because it does not favor more complete structures. In addition, it is
        very sensitive to minor changes in number of basepairs and nts.

        Parameters
        ----------
        members : dict
            A group dict to find the representative of.

        Returns
        --------
            The initial representative.
        """
        naive = Naive(self.config, self.session)
        return naive(group)[0]

    def candidates(self, best, members):
        """Find all possible candidates for a representative within the group,
        given a current best ife. This finds all chains that have at least as
        many basepairs and nucleotides as the best chain. The chains will be
        returned in sorted order.

        :param dict best: The current best.
        :param list group: The list of dicts to search.
        :returns: The list of candidates for the representative.
        """

        len = op.itemgetter('length')
        bp = op.itemgetter('bp')
        same = lambda c: bp(c) == bp(best) and len(c) == len(best)
        possible = it.ifilter(lambda c: len(c) >= len(best), members)
        possible = it.ifilter(lambda c: bp(c) >= bp(best), possible)
        possible = it.ifilterfalse(same, possible)
        return sorted(possible, key=bp_per_nt, reverse=True)

    def increase(self, first, second, key):
        """Compute the percent increase for the given set of dictionaries and
        with the given key. If the second one is 0 then we return 100 for 100%
        increase.

        :param dict first: Dictionary to get the increase to.
        :param dict second: Dictionary to get the increase from.
        :param str key: Key to use
        :returns: The percent increase.
        """

        if not second[key]:
            if not first[key]:
                return 0
            return 100
        return (float(first[key]) / float(second[key]) - 1) * 100

    def best_above_cutoffs(self, representative, candidates, length_increase,
                           bp_increase):
        """This will find the true representative given a current one and a
        list of candidates. This will attempt to maximize the number of bps and
        nts in the candidate as compared to the current representative. In
        addition, it will only change representatives if we have have enough of
        an increase. This adds stability to the process so minor improvements
        are ignored, while large ones will lead to large changes.

        :param dict representative: The current representative.
        :param list candidates: A list of candidates to examine.
        :returns: The new representative.
        """

        for candidate in candidates:
            length_change = self.increase(candidate, representative, 'length')
            bp_change = self.increase(candidate, representative, 'bp')
            if length_change >= length_increase and bp_change >= bp_increase:
                return candidate
        return representative

    def __call__(self, initial, length_increase=NR_LENGTH_PERCENT_INCREASE,
                 bp_increase=NR_BP_PERCENT_INCREASE):
        """Find the representative for the group.

        Parameters
        ----------
        group : list
            List of IFE's to find the representative of.
        length_increase : float
            The fraction increase in resolved nucleotides that an IFE must have
            to be selected as representative.
        bp_increase : float
            The fraction increase of basepairs that an IFE must have to be
            selected as representative.

        Returns
        -------
            representative : dict
        The ife which should be the representative.
        """

        group = self.filter_group_by_method(initial)
        best = self.initial_representative(group)
        if not best:
            raise core.InvalidState("No current representative")
        self.logger.debug("Naive representative: %s", best['id'])

        rep = best
        while True:
            candidates = self.candidates(rep, group['members'])
            self.logger.debug("Found %i representative candidates",
                              len(candidates))
            new_rep = self.best_above_cutoffs(rep, candidates,
                                              length_increase,
                                              bp_increase)
            if new_rep == rep:
                break
            self.logger.info("Changed representative from %s to %s",
                             rep['id'], new_rep['id'])
            rep = new_rep

        if not rep:
            raise core.InvalidState("No representative found")

        self.logger.debug("Computed representative: %s", rep['id'])
        return self.insert_as_representative(rep,
                                             initial['members'],
                                             sort=bp_per_nt)


class AnyIncrease(Increase):
    """A modification of the percent-increase representative method that will
    instead select the representative if it has any increase in length and bp.
    """
    method = 'any-increase'

    def __call__(self, group):
        return super(AnyIncrease, self).__call__(group,
                                                 length_increase=0,
                                                 bp_increase=0)


class ParentIncrease(Increase):
    """This is a modification of the percent increase method which uses the
    parent representative as the initial representative, if possible. If this
    is not possible then it will go back to the simple bp/nt selection.
    """
    method = 'parent-percent-increase'

    def initial_representative(self, group):
        """Select the initial representative. This will be the representative
        for the parent if there is only 1 parent. If there is more than 1
        parent, or none, then this will switch to bp/nt.

        Parameters
        ----------
        group : dict
            A group dict that must have a parent entry.

        Returns
        -------
            The initial representative.
        """

        if not group['parent']:
            self.logger.debug("No parent to use representative of")
            return super(ParentIncrease, self).initial_representative(group)

        if len(group['parent']) > 1:
            self.logger.debug("Multiple parents, using naive")
            return super(ParentIncrease, self).initial_representative(group)

        parent = group['parent'][0]
        previous = parent.get('representative', None)
        if not previous:
            self.logger.warning("Loaded parent %s has no representative, using"
                                " naive", parent)
            return super(ParentIncrease, self).initial_representative(group)

        current_ids = {mem['id'] for mem in group['members']}
        if previous['id'] not in current_ids:
            self.logger.warning("Not using previous rep %s, because not in "
                                "current group: %s",
                                previous['id'], current_ids)
            return super(ParentIncrease, self).initial_representative(group)

        return previous


class QualityMetrics(Representative):
    """Find representatives using quality metrics provided by PDB.
    """
    method = 'quality-metrics'
    hardcoded = MANUAL_IFE_REPRESENTATIVES
    worse = WORSE_THAN_MANUAL_IFE_REPRESENTATIVES

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
                                  mod.PdbQuality.percent_rsrz_outliers.label('rsrz'),
                                  mod.PdbQuality.clashscore,
                                  mod.PdbQuality.percent_rota_outliers.label('backbone'),
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
            required = set(['rsrz', 'backbone', 'clashscore'])
            return member['method'] == 'X-RAY DIFFRACTION' and \
                member['resolution'] <= 4.0 and \
                required.issubset(member['quality']['has'])

        if any(is_good_xray(m) for m in members):
            return [m for m in members if is_good_xray(m)]
        return list(members)

    def filter_by_nts(self, members):
        best = max(m['length'] for m in members)
        return [m for m in members if m['length'] >= 0.75 * best]

    def filter_by_resolution(self, members):
        best = min(m['resolution'] for m in members)
        return [m for m in members if abs(m['resolution'] - best) <= 0.2]

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
        best_method = self.filter_by_method(with_quailty)
        best_nts = self.filter_by_nts(best_method)
        best_resolution = self.filter_by_resolution(best_nts)
        if not best_resolution:
            raise core.InvalidState("Nothing with good resolutin")
        ordered_by_quality = self.sort_by_quality(best_resolution)
        with_representative = self.use_hardcoded(ordered_by_quality)
        return self.final_ordering(with_representative, given['members'])
