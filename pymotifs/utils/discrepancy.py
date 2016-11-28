"""This contains some utility functions for dealing with discrepancies.
"""

from pymotifs.constants import MAX_RESOLUTION_DISCREPANCY
from pymotifs.constants import MIN_NT_DISCREPANCY


def should_compare_chain_discrepancy(chain):
    """Check if we can compared discrepancies using this chain.

    Parameters
    ----------
    chain : dict
        The chain dict to test.

    Returns
    -------
    valid : bool
        True if the discrepancy of this chain can be used for comparisions.
    """
    return valid_chain(chain)


def should_compute_chain_discrepancy(chain):
    """Check if we should compute the discrepancy using this chain.

    Parameters
    ----------
    chain : dict
        The chain dict to test.

    Returns
    -------
    valid : bool
        True if this chain should have a discrepancy computed using it.
    """
    return valid_chain(chain)


def valid_chain(chain):
    """Check if the chain can have a dsicrepancy computed. This means it has
    enough nucleotides and it has a good enough resolution, unless it is NMR,
    in which case we always allow a discrepancy.

    Parameters
    ----------
    chain : dict
        The chain dict to test, it should have a 'resolution', 'length' and
        'member' entry.

    Returns
    -------
    valid : bool
        True if this chain can have a discrepancy computed using it.
    """

    if 'length' not in chain:
        return False

    if chain['length'] < MIN_NT_DISCREPANCY:
        return False

    if 'method' not in chain:
        return False

    if chain['method'] != 'SOLUTION NMR':
        return chain['resolution'] is not None and \
            chain['resolution'] <= MAX_RESOLUTION_DISCREPANCY
    return True
