import collections as coll

from pymotifs import models as mod
from pymotifs import core
from pymotifs.utils import structures as st


def autonomous_sorter(chain):
    # When autonomous is set to None, this will default to False
    autonomous = chain.get('autonomous') or False
    return (int(autonomous), chain['length'], -ord(chain['name']))


class Info(core.Base):
    def __init__(self, *args, **kwargs):
        super(Info, self).__init__(*args, **kwargs)
        self.bp_helper = st.BasePairQueries(self.session.maker)
        self.struct_helper = st.Structure(self.session.maker)

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

    def load(self, pdb, chain, molecule_type='rna', use_names=False):
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
            'autonomous': None
        }

        with self.session() as session:
            query = session.query(mod.ChainInfo.chain_id,
                                  mod.ChainInfo.chain_length,
                                  mod.ChainInfo.entity_name,
                                  mod.ChainInfo.sequence).\
                filter_by(pdb_id=pdb, chain_name=chain)

            result = query.one()
            data['db_id'] = result.id
            data['exp_length'] = result.chain_length
            data['entity'] = result.entity_name
            data['sequence'] = result.sequence

        data['source'] = None
        try:
            data['source'] = self.struct_helper.source(pdb, chain,
                                                       simplify=True)
        except:
            self.logger.warn("Failed to find all taxon ids for %s", data['id'])

        with self.session() as session:
            query = session.query(mod.UnitInfo.unit_id).\
                filter_by(pdb_id=pdb, chain=chain, unit_type_id=molecule_type)
            data['length'] = query.count()

        data.update(self.bps(pdb, chain))
        return data

    def load_all(self, pdb, **kwargs):
        """Load all RNA containing chains in a pdb file.

        :pdb: The PDB id.
        """

        rna_chains = self.rna_chains(pdb)
        return [self.load(pdb, chain, **kwargs) for chain in rna_chains]

    def cross_chain_interactions(self, chains):
        """Create a dictionary of the interactions between the listed chains.
        This will get only the counts.

        :chains: A list of chain dictionaries.
        :returns: A dictionary of like { 'A': { 'B': 10 }, 'B': { 'A': 10 } }.
        """

        helper = st.BasePairQueries(self.session.maker)
        pdb = chains[0]['pdb']
        interactions = coll.defaultdict(dict)

        for chain1 in chains:
            name1 = chain1['name']
            for chain2 in chains:
                if chain1 == chain2:
                    continue

                name2 = chain2['name']
                if name2 in interactions.get(name1, {}):
                    continue

                count = helper.cross_chain(pdb, name1, name2, count=True)
                interactions[name1][name2] = count
                interactions[name2][name1] = count

        return interactions

    def merge(self, chains):
        """This merges a list of chain dictionaries into one dictionary. It will
        contain the same information as the original chains but merged for all
        chains. This means that the id is a comma separated string of the
        composite ids. The pdb is still a string. db_id, name, and entity are
        lists of all the chain information. Internal, bp and length are the sum
        of the given list. External is recomputed for this set of chains.

        :chains: A list of chains to merge.
        :returns: The merged dictionary.
        """

        helper = self.bp_helper
        merged = {
            'id': [],
            'pdb': chains[0]['pdb'],
            'internal': 0,
            'external': None,
            'bp': 0,
            'length': 0,
            'exp_length': 0,
            'lr': 0,
            'count': len(chains)
        }

        sorted_chains = list(reversed(sorted(chains, key=autonomous_sorter)))
        for chain in sorted_chains:
            merged['id'].append(chain['id'])
            merged['internal'] += chain['internal']
            merged['bp'] += chain['bp']
            merged['length'] += chain['length']
            merged['exp_length'] += chain['exp_length']
            merged['lr'] += chain['lr']

        merged['id'] = ','.join(merged['id'])
        names = [chain['name'] for chain in chains]
        merged['external'] = helper.cross_chain(merged['pdb'], names,
                                                count=True, family='cWW')

        return {
            'id': merged.pop('id'),
            'pdb': merged.pop('pdb'),
            'summary': merged,
            'chains': sorted_chains
        }
