import inspect
import operator as op
import itertools as it

from pymotifs import core

from pymotifs.constants import NR_BP_PERCENT_INCREASE
from pymotifs.constants import NR_LENGTH_PERCENT_INCREASE


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
    pass


def known():
    finders = []
    for key, value in globals().items():
        if inspect.isclass(value) and issubclass(value, Representative) and \
                value is not Representative:
            finders.append((value.name, value))
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
                value.name == name:
            return value
    raise ValueError("Unknown method %s" % name)


class Naive(Representative):
    """This will select the naive representative. This just selects the chain
    with the most bp/nt and uses that as the representative.
    """
    name = 'naive'

    def __call__(self, group):
        return sorted(group, key=bp_per_nt, reverse=True)[0]


class Increase(Representative):
    """A class to find the representative for a group of ifes. This will find
    the best in terms of bp/nt and attempt to find all those with more bp's and
    nts in the set.
    """
    name = 'increase'

    def naive_best(self, group):
        """Find the best chain terms of bps/nts. This is the starting point for
        finding the representative in a set of ifes. This method is naive
        because it does not favor more complete structures. In addition, it is
        very sensitive to minor changes in number of basepairs and nts.

        :param list group: A list of dictonaries to find the naive
        representative of.
        :returns: The initial representative.
        """
        return max(group, key=bp_per_nt)

    def candidates(self, best, group):
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
        possible = it.ifilter(lambda c: len(c) >= len(best), group)
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

    def best_above_cutoffs(self, representative, candidates):
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

        cutoff = (NR_LENGTH_PERCENT_INCREASE, NR_BP_PERCENT_INCREASE)
        possible = []
        for candidate in candidates:
            length_change = self.increase(candidate, representative, 'length')
            bp_change = self.increase(candidate, representative, 'bp')
            if (length_change, bp_change) >= cutoff:
                possible.append(candidate)

        if not possible:
            return representative
        return max(possible, key=bp_per_nt)

    def __call__(self, possible):
        """Find the representative for the group.

        :param list group: List of ifes to find the best for.
        :returns: The ife which should be the representative.
        """

        # Prefer any xray over any cyro em, as cyro modesl are generally built
        # using x-ray and not yet carefully modeled.
        group = [ife for ife in possible if ife['method'] == 'X-RAY DIFFRACTION']
        if not group:
            group = possible

        best = self.naive_best(group)
        if not best:
            raise core.InvalidState("No current representative")
        self.logger.debug("Naive representative: %s", best['id'])

        candidates = self.candidates(best, group)
        self.logger.debug("Found %i representative candidates",
                          len(candidates))

        rep = self.best_above_cutoffs(best, candidates)
        if not rep:
            raise core.InvalidState("No representative found")

        if rep['id'] != best['id']:
            self.logger.info("Changed representative from %s to %s",
                             best['id'], rep['id'])

        self.logger.debug("Computed representative: %s", rep['id'])
        return rep
