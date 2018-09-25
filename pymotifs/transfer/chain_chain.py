"""A module to handle importing and exporting chain chain data. Sometimes we
want to copy chain chain data from one machine to another. This module is meant
to have the logic for doing so in one place.
"""

import pickle
import logging
import operator as op

from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import aliased

from pymotifs.core import Session
from pymotifs import models as mod
from pymotifs.utils import row2dict

chain_entry = lambda n: op.itemgetter('pdb_id' + n, 'chain_name' + n)
chain1 = chain_entry('1')
chain2 = chain_entry('2')

logger = logging.getLogger(__name__)


def setup(engine, **kwargs):
    """Reflect the models and create a session wrapper.

    Parameters
    ----------
    engine : Engine
        An engine to use to reflect models.

    Returns
    session : pymotifs.core.Session
        The session wrapper to use.
    """
    mod.reflect(engine)
    return Session(sessionmaker(engine))


def chain_id_mapping(session, data, ignore_missing=False):
    """Create a mapping from (pdb, chain_name) to chain id in the database.
    By default, this will map all entries in data or it will fail with a value
    error.

    Parameters
    ----------
    session : pymotifs.core.Session
        The session to use
    data : list
        A list of dictionaries with pdb_id1, chain_name1, pdb_id2, chain_name2
        keys.

    Returns
    -------
    mapping : dict
        A dict mapping from (pdb, chain_name) to chain id for all entries in
        data.
    """

    mapping = {}
    entries = set()
    entries.update(chain1(e) for e in data)
    entries.update(chain2(e) for e in data)
    with session() as sess:
        query = sess.query(mod.ChainInfo)
        for result in query:
            entry = (result.pdb_id, result.chain_name)
            if entry not in entries:
                continue
            mapping[entry] = result.chain_id
            entries.remove(entry)

    if entries:
        msg = "Not all entries mapped to chain ids: %s" % str(entries)
        logger.error(msg)
        if not ignore_missing:
            raise ValueError(msg)

    return mapping


def correspondence_id_mapping(session, data, ignore_missing=False):
    """Create a mapping from compared chain chain to correspondence ids. This
    will fail if not all chains in the input data could be mapped, if
    ignore_missing is False (the default behavior), otherwise it will only log
    the error.

    Parameters
    ----------
    session : pymotifs.core.Session
        The sesson to use
    data : list
        A list of dictionaries with pdb_id1, chain_name1, pdb_id2, chain_name2
        entries.
    ignore_missing : bool, optional
        A flag to make this ignore missing chains. In this case errors are only
        logged.

    Returns
    -------
    mapping : dict
        A dictionary mapping (chain_id, chain_id) to correspondence id.
    """

    entries = {(chain1(e), chain2(e)) for e in data}

    with session() as sess:
        corr = mod.CorrespondencePdbs
        query = sess.query(corr.correspondence_id,
                           corr.pdb_id_1.label('pdb_id1'),
                           corr.pdb_id_2.label('pdb_id2'),
                           corr.chain_name_1.label('chain_name1'),
                           corr.chain_name_2.label('chain_name2'),
                           corr.chain_id_1,
                           corr.chain_id_2,
                           )

        mapping = {}
        for result in query:
            result = row2dict(result)
            ids = (chain1(result), chain2(result))
            if ids not in entries:
                continue
            key = (result['chain_id_1'], result['chain_id_2'])
            if key in mapping:
                raise ValueError("Duplicate mapping found %s" % ids)
            mapping[key] = result['correspondence_id']
            entries.remove(ids)

    logger.info("Found %i/%i correspondences", len(mapping), len(entries))
    if entries:
        logger.error("Could not map all correspondences %s", str(entries))
        if not ignore_missing:
            raise ValueError("Could not map all correspondences %s" %
                             str(entries))

    return mapping


def known_comparisions(session):
    """Get a set of all known chain chain comparisons. This can be used to skip
    importing existing comparisons.

    Parameters
    ----------
    session : pymotifs.core.Session
        The session to use.

    Returns
    -------
    known : set
        A set of all known (chain_id, chain_id) comparisons.
    """

    known = set()
    with session() as sess:
        query = sess.query(mod.ChainChainSimilarity)
        for result in query:
            known.add((result.chain_id_1, result.chain_id_2))
    return known


def load(filename, ignore_missing=False, **kwargs):
    """Load the data in the given file. This will import all new data in the
    given file to the database. If the comparison has already been done it will
    not be loaded. The file should contain the data to write in the format
    produced by dump.

    Parameters
    ----------
    filename : str
        The name of the file to load data from.

    Raises
    ------
    ValueError
        If either any chain is not mapped or there is not correspondence
        between a compared pair of chains and ignore_mapping is False.
    """

    with open(filename, 'rb') as raw:
        data = pickle.load(raw)

    session = setup(**kwargs)
    chain_mapping = chain_id_mapping(session, data,
                                     ignore_missing=ignore_missing)
    corr_mapping = correspondence_id_mapping(session, data,
                                             ignore_missing=ignore_missing)
    known = known_comparisions(session)
    savable = []
    for entry in data:
        if chain1(entry) not in chain_mapping:
            logger.error("Chain %s not in chain mapping", chain1(entry))
            if not ignore_missing:
                raise ValueError("Missing chain id")

        if chain2(entry) not in chain_mapping:
            logger.error("Chain %s not in chain mapping", chain2(entry))
            if not ignore_missing:
                raise ValueError("Missing chain id")

        chain_id_1 = chain_mapping[chain1(entry)]
        chain_id_2 = chain_mapping[chain2(entry)]
        ids = (chain_id_1, chain_id_2)
        if ids in known:
            logger.info("Not importing known correspondence between %s %s",
                        chain1(entry), chain2(entry))
            continue

        if ids not in corr_mapping:
            logger.error("Missing correspondence between %s", ids)
            if not ignore_missing:
                raise ValueError("Missing correspondence")
            continue

        to_save = dict(entry)
        del to_save['pdb_id1']
        del to_save['pdb_id2']
        del to_save['chain_name1']
        del to_save['chain_name2']
        to_save.update({
            'chain_id_1': chain_id_1,
            'chain_id_2': chain_id_2,
            'correspondence_id': corr_mapping[ids],
        })
        logger.debug("Importing %s", to_save)
        savable.append(mod.ChainChainSimilarity(**to_save))

    with session() as sess:
        for entry in savable:
            sess.insert(entry)
        logger.info("Imported %i correspondences", len(savable))


def dump(filename, **kwargs):
    """Dump chain chain comparison data to a file. This will dump all chain
    chain comparison data to a file for later import. The data is pickled for
    easy reading and writing in python.

    Parameters
    ----------
    filename : str
        Name of the file to write to.
    """

    session = setup(**kwargs)
    with session() as sess:
        chain1 = aliased(mod.ChainInfo)
        chain2 = aliased(mod.ChainInfo)
        query = sess.query(mod.ChainChainSimilarity.discrepancy,
                           mod.ChainChainSimilarity.num_nucleotides,
                           mod.ChainChainSimilarity.model_1,
                           mod.ChainChainSimilarity.model_2,
                           chain1.pdb_id.label('pdb_id1'),
                           chain1.chain_name.label('chain_name1'),
                           chain2.pdb_id.label('pdb_id2'),
                           chain2.chain_name.label('chain_name2'),
                           ).\
            join(chain1,
                 chain1.chain_id == mod.ChainChainSimilarity.chain_id_1).\
            join(chain2,
                 chain2.chain_id == mod.ChainChainSimilarity.chain_id_2)

        results = [row2dict(r) for r in query]
        with open(filename, 'wb') as out:
            pickle.dump(results, out)
