import itertools as it
import collections as coll

from pymotifs import core
from pymotifs import models as mod
from pymotifs.nr import connectedsets as cs
from pymotifs.utils import structures as st
from pymotifs.utils import correspondence as cr


class Grouper(core.Stage):
    """This is a class to determine the grouping of chains in a pdb. This will
    load all non-redundant chains and then group them by similarity. It
    produces a list of chains which are non-redundant. This does not process
    all chains as it only takes the best of each type from files. However, it
    does make a note about those other chains.
    """

    internal_cutoff = 0.5

    def __init__(self, *args, **kwargs):
        super(Grouper, self).__init__(*args, **kwargs)
        self.bp_helper = st.BasePairQueries(self.session.maker)
        self.struct_helper = st.Structure(self.session.maker)

    def best(self, chains):
        """Compute the best chain for a list of chains

        :chains: A list of chains to select the best one of.
        :returns: A single entry from the list that is the best. It adds an
        entry to the chains called 'equivalent', which is all other chains in
        that structure which are equivalent.
        """

        def ordering(chain):
            if chain['bp'] == 0:
                bp_nt = 0
            else:
                bp_nt = float(chain['bp']) / float(chain['length'])

            name = chain['name']
            if isinstance(name, list):
                name = min(name)

            return (bp_nt, chain['length'], name)

        best = max(chains, key=ordering)
        best['equivalent'] = list(chains)
        best['equivalent'].remove(best)

        return best

    def rna_chains(self, pdb):
        """Get all RNA containing chains for a given pdb.

        :pdb: The PDB file to get chains for.
        :returns: A list of chain ids for that pdb.
        """
        with self.session() as session:
            query = session.query(mod.ChainInfo.chain_name).\
                filter_by(entity_macromolecule_type='Polyribonucleotide (RNA)').\
                filter_by(pdb_id=pdb)
            return [result.chain_name for result in query]

    def merge_chains(self, chains):
        """This merges a list of chain dictonaries into one dictonary. It will
        contain the same information as the original chains but merged for all
        chains. This means that the id is a comma seperated string of the
        composite ids. The pdb is still a string. db_id, name, and entity are
        lists of all the chain information. internal, bp and length are the sum
        of the given list. external is recomputed for this set of chains.

        :chains: A list of chains to merge.
        :returns: The merged dictonary.
        """

        helper = self.bp_helper
        merged = {
            'id': [],
            'pdb': chains[0]['pdb'],
            'db_id': [],
            'name': [],
            'internal': 0,
            'external': None,
            'bp': 0,
            'length': 0,
            'exp_length': 0,
            'entity': [],
            'source': [],
            'lr': 0
        }
        chains.sort(key=lambda c: c['name'])

        for chain in chains:
            merged['id'].append(chain['id'])
            merged['db_id'].append(chain['db_id'])
            merged['name'].append(chain['name'])
            merged['entity'].append(chain['entity'])
            merged['source'].append(chain['source'])
            merged['internal'] += chain['internal']
            merged['bp'] += chain['bp']
            merged['length'] += chain['length']
            merged['exp_length'] += chain['exp_length']
            merged['lr'] += chain['lr']

        merged['id'] = ','.join(merged['id'])
        merged['entity'].sort()
        merged['external'] = helper.cross_chain(merged['pdb'], merged['name'],
                                                count=True, family='cWW')

        if set(merged.keys()) != set(chains[0].keys()):
            missing = set(chains[0].keys()) - set(merged.keys())
            self.logger.warning("Did not merge with all keys. Missing: %s",
                                missing)

        return merged

    def autonomous(self, chains):
        """This method goes through all chains and detects those that are
        autonomous. When a chain is making no interactions with other chains it
        is kept, if the chain is making important interactions with another
        chain then those will be taken as one autonomous unit and kept.

        :chains: A list of chain dictionaries.
        :returns: A list of merged chain objects.
        """

        if len(chains) == 1:
            return [self.merge_chains(chains)]

        mapping = {}
        grouped = []
        connections = coll.defaultdict(set)
        helper = self.bp_helper

        for chain in chains:
            pdb = chain['pdb']
            chain_id = chain['name']
            mapping[chain['id']] = chain

            if chain['internal']:
                internal = float(chain['internal'])
                external = float(chain['external'])
                total = internal + external
                if internal / total >= self.internal_cutoff:
                    grouped.append(self.merge_chains([chain]))
                    continue

            for chain2 in chains:
                if chain2 == chain:
                    next
                if helper.cross_chain(pdb, chain_id, chain2['name']):
                    connections[chain['id']].add(chain2['id'])

            if chain['id'] not in connections:
                connections[chain['id']] = set()

        groups = cs.find_connected(dict(connections))

        for group in groups.values():
            merged = self.merge_chains([mapping[name] for name in group])
            grouped.append(merged)

        return grouped

    def bps(self, pdb, chain):
        """Determine the total basepairs, internal cWW, external cWW,
        and long range basepairs.

        :pdb: A pdb id to lookup.
        :chain: A chain to use.
        :returns: A dictionary with bp, internal, external and lr entries,
        which stores the respective data.
        """

        helper = self.bp_helper
        return {
            'bp': helper.representative(pdb, chain, count=True),
            'lr': helper.representative(pdb, chain, count=True,
                                        range_cutoff=st.LONG_RANGE_CUTOFF),
            'internal': helper.representative(pdb, chain, count=True,
                                              family='cWW'),
            'external': helper.cross_chain(pdb, chain, count=True,
                                           family='cWW')
        }

    def info(self, pdb, chain, molecule_type='rna', use_names=False):
        """This loads all information about a chain into a dictionary. This
        will load generic information about a chain, such as resolved, length,
        database id, the source and information about basepairing. The
        basepairing information is loaded from `bps`.

        :pdb: The pdb to search.
        :chain: The chain to search.
        :returns: A dictionary with
        """

        data = {
            'id': '%s|%s' % (pdb, chain),
            'name': chain,
            'pdb': pdb,
        }

        with self.session() as session:
            query = session.query(mod.ChainInfo.id,
                                  mod.ChainInfo.chain_length,
                                  mod.ChainInfo.entity_name).\
                filter_by(pdb_id=pdb, chain_name=chain)

            result = query.one()
            data['db_id'] = result.id
            data['exp_length'] = result.chain_length
            data['entity'] = result.entity_name

        data['source'] = None
        try:
            data['source'] = self.struct_helper.source(pdb, chain,
                                                       simplify=True)
        except:
            self.logger.warn("Failed to find all taxon ids for %s", data['id'])

        with self.session() as session:
            query = session.query(mod.UnitInfo.id).\
                filter_by(pdb_id=pdb, chain=chain, unit_type_id=molecule_type)
            data['length'] = query.count()

        data.update(self.bps(pdb, chain))
        return data

    def chains(self, pdb):
        """Load all chains from a given pdb. This will get the RNA chains as
        well as load some interaction data about the chains.

        :pdb: A pdb to get the chains for.
        :returns: A list of dictionaries with data about all chains. The data
        is that which is provided by the info method.
        """

        chains = [self.info(pdb, chain) for chain in self.rna_chains(pdb)]
        return self.autonomous(chains)

    def nr_chains(self, pdb):
        """Get all chains which are non-redundant from a cif file.
        """
        chains = self.chains(pdb)
        key = lambda c: c['entity']
        grouped = it.groupby(sorted(chains, key=key), key)
        best = []
        for _, group in grouped:
            best.append(self.best(list(group)))
        return best

    def alignment_info(self, chains):
        """Load the data about alignments between all pairs of chains.

        :chains: A list of chains.
        :returns: A dictionary of dictionaries.
        """
        pdbs = set(chain['pdb'] for chain in chains)

        helper = cr.Helper(self.session.maker)
        return helper.aligned_chains(pdbs)

    def is_equivalent(self, alignments, chain1, chain2):
        """Determine if two chains are equivalent.

        :alignments: A dictionary of dictionaries about the alignment between
        chains.
        :chain1: The first chain.
        :chain2: The second chain.
        :returns: True if the two chains are equivalent, False otherwise
        """
        if len(chain1['source']) != len(chain2['source']):
            return False

        for source1, source2 in zip(chain1['source'], chain2['source']):
            if source1 and source2 and source1 != source2:
                self.logger.debug("Splitting %s, %s: Conflicting species",
                                  chain1['id'], chain2['id'])
                return False

        # TODO: Split only if entire alignment is bad. This has a different set
        # of cutoffs as the a small chain may not align but a large may leading
        # to an overall valid alignment
        ids1 = it.chain(chain1['db_id'])
        ids2 = it.chain(chain2['db_id'])
        for id1, id2 in zip(ids1, ids2):
            if id1 not in alignments:
                self.logger.debug("Splitting %s %s: Unaligned %s",
                                  chain1['id'], chain2['id'], chain1['id'])
                return False

            if id2 not in alignments[id1]:
                self.logger.debug("Splitting %s, %s: Poor alignment",
                                  chain1['id'], chain2['id'])
                return False

        return True

    def group(self, chains, alignments):
        """Group all chains into connected components.

        :chains: A list of chains to group.
        :returns: A list of lists of the connected components.
        """

        graph = coll.defaultdict(set)
        mapping = {}
        for chain1 in chains:
            mapping[chain1['id']] = chain1

            for chain2 in chains:
                if chain1 == chain2 or \
                        self.is_equivalent(alignments, chain1, chain2):
                    graph[chain1['id']].add(chain2['id'])
                    graph[chain2['id']].add(chain1['id'])

        groups = []
        for ids in cs.find_connected(graph).values():
            groups.append(mapping[id] for id in ids)
        return groups

    def representative(self, group):
        """Compute the representative for a group of chains.

        :group: A list of chains.
        :returns: The representative entry from the list.
        """
        rep = self.best(group)
        del rep['equivalent']
        return rep

    def __call__(self, pdbs, **kwargs):
        chains = it.imap(self.nr_chains, pdbs)
        chains = list(it.chain.from_iterable(chains))

        alignments = self.alignment_info(chains)
        groups = []
        for group in self.group(chains, alignments):
            members = list(group)
            representative = self.representative(members)
            members.remove(representative)
            groups.append({
                'members': members,
                'representative': representative
            })

        return groups
