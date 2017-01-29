import abc
import copy
import inspect
import operator as op
import itertools as it

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import known_subclasses
from pymotifs.utils import row2dict

from pymotifs.constants import NR_BP_PERCENT_INCREASE
from pymotifs.constants import NR_LENGTH_PERCENT_INCREASE
from pymotifs.constants import NR_ALLOWED_METHODS
from pymotifs.constants import NR_HARD_CODED


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
        return sorted(group['members'], key=bp_per_nt, reverse=True)[0]


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
        return naive(group)

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
        return rep


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

    def hardcoded_representative(self, group):
        members = {m['id'] for m in group['members']}
        found = members.intersection(NR_HARD_CODED)
        if len(found) > 1:
            self.logger.error("More than one hardcoded, using quality for %s",
                              group)
            return None

        if not found:
            self.logger.debug("No hardcoded representative for %s", group)
            return None
        return found.pop()

    def select_final_representative(self, group):
        selector = ParentIncrease(self.config, self.session)
        return selector(group)

    def lookup_structure_quality(self, members):
        pdbs = set(member['pdb'] for member in members)
        with self.session() as session:
            query = session(mod.PdbQuality).\
                filter(mod.PdbQuality.pdb_id.in_(pdbs))
            measures = {}
            for result in query:
                measures[result.pdb_id] = row2dict(result)
        return measures

    def lookup_chain_quality(self, members):
        pdbs = set(member['pdb'] for member in members)
        with self.session() as session:
            query = session(mod.ChainQuality).\
                filter(mod.ChainQuality.pdb_id.in_(pdbs))
            measures = {}
            for result in query:
                measures[result.chain_id] = row2dict(result)
        return measures

    def filter_by_quality(self, group):
        return group

    def __call__(self, given):
        hardcoded = self.hardcoded_representative(given)
        if hardcoded:
            return hardcoded

        if self.is_cryo_group(given):
            return self.filter_by_bp(given)

        group = self.filter_group_by_method(given)
        group = self.filter_group_by_quality(group)
        return self.select_final_representative(group)
