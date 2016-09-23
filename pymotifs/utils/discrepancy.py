"""This contains some utility functions for dealing with discrepancies.
"""

from pymotifs.constants import MAX_RESOLUTION_DISCREPANCY
from pymotifs.constants import MIN_NT_DISCREPANCY


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

    if chain['length'] < MIN_NT_DISCREPANCY:
        return False

    if chain['method'] != 'SOLUTION NMR':
        return chain['resolution'] is not None and \
            chain['resolution'] <= MAX_RESOLUTION_DISCREPANCY
    return True
