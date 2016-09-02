"""A module to handle importing and exporting chain chain data. Sometimes we
want to copy chain chain data from one machine to another. This module is meant
to have the logic for doing so in one place.
"""

import sys
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
    """Create a session wrapper.
    """
    mod.reflect(engine)
    return Session(sessionmaker(engine))


def chain_id_mapping(session, data, ignore_missing=False):
    """Create a mapping from (pdb, chain_name) to chain id in the database.
    By default, this will map all entries in data or it will fail with a value
    error.

    :param Session session: The session to use.
    :param list data: A list of dictionaries with pdb_id1, chain_name1,
    pdb_id2, chain_name2 keys.
    :returns: A dict mapping from (pdb, chain_name) to chain id for all entries
    in data.
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

    if entries and not ignore_missing:
        raise ValueError("Not all entries mapped: %s" % str(entries))

    return mapping


def correspondence_id_mapping(session, data, ignore_missing=False):
    """Create a mapping from compared chain chain to correspondence ids. This
    will fail if not all chains in the input data could be mapped.

    :param Session session: The session to use.
    :param list data: A list of dictionaries with pdb_id1, chain_name1, pdb_id2, chain_name2
    entries.
    :param Bool ignore_missing: A flag to make this ignore missing chains. In
    this case errors are logged instead.
    :returns: A dictionary mapping (chain_id, chain_id) to correspondence id.
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
        self.logger.error("Could not map all correspondences %s", str(entries))
        if not ignore_missing:
            raise ValueError("Could not map all correspondences %s" % str(entries))

    return mapping


def known_comparisions(session):
    """Get a set of all known chain chain comparisons. This can be used to
    skip importing existing comparisons.

    :param Session session: The session to use.
    :returns: A set of all known (chain_id, chain_id) comparisons.
    """

    known = set()
    with session() as sess:
        query = sess.query(mod.ChainChainSimilarity)
        for result in query:
            known.add((result.chain_id_1, result.chain_id_2))
    return known


def load(filename, **kwargs):
    """Load the data in the given file. This will import all new data in the
    given file to the database. If the comparison has already been done it
    will not be loaded. The file should contain the data to write in the format
    produced by dump.

    :param str filename: The name of the file to load data from.
    """

    with open(filename, 'rb') as raw:
        data = pickle.load(raw)

    session = setup(**kwargs)
    chain_mapping = chain_id_mapping(session, data)
    correspondence_mapping = correspondence_id_mapping(session, data)
    known = known_comparisions(session)
    savable = []
    for entry in data:
        chain_id_1 = chain_mapping[chain1(entry)]
        chain_id_2 = chain_mapping[chain2(entry)]
        ids = (chain_id_1, chain_id_2)
        if ids in known:
            logger.info("Not importing known correspondence between %s %s",
                        chain1(entry), chain2(entry))
            continue

        to_save = dict(entry)
        del to_save['pdb_id1']
        del to_save['pdb_id2']
        del to_save['chain_name1']
        del to_save['chain_name2']
        to_save.update({
            'chain_id_1': chain_id_1,
            'chain_id_2': chain_id_2,
            'correspondence_id': correspondence_mapping[ids],
        })
        logger.debug("Importing %s", to_save)
        savable.append(mod.ChainChainSimilarity(**to_save))

    with session() as sess:
        for entry in savable:
            sess.insert(entry)
        logger.info("Imported %i correspondences", len(savable))


def dump(filename, **kwargs):
    """Dump chain chain comparison data to a file. This will dump all chain
    chain comparison data to a file for later import.

    :param str filename: Name of the file to write to.
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
