"""Map chains to experimental sequences. This will only process RNA chains and
will not map protein chains to to experimental sequences as we only compute
experimental sequence data for RNA chains.
"""

import operator as op
import itertools as it
import functools as ft

from pymotifs import core

from pymotifs import models as mod

from pymotifs.utils.structures import Structure

from pymotifs.chains.info import Loader as ChainLoader
from pymotifs.exp_seq.info import Loader as InfoLoader

from sqlalchemy import case


class Loader(core.SimpleLoader):
    dependencies = set([ChainLoader, InfoLoader])    

    def to_process_old(self, pdbs, **kwargs):
        """Compute all chain ids to process. This will extract all rna chain
        ids in the given list of pdbs.

        Parameters
        ----------
        pdbs : list
            A list of `str` pdb ids to process

        Returns
        -------
        chain_ids : list
            A list of sorted `int` chain ids.
        """

        helper = Structure(self.session.maker)
        na_chains = ft.partial(helper.na_chains, return_id=True)
        chains = it.imap(na_chains, pdbs)
        chains = it.chain.from_iterable(chains)
        chains = it.imap(op.itemgetter(1), chains)

        return sorted(set(chains))

    def to_process(self,pdbs,**kwargs):
        """Compute current chain ids to process. This will extract current chain
        ids in the given list of pdbs.
        Parameters
        ----------
        pdbs : list
            A list of `str` pdb ids to process
        Returns
        -------
        chain_id : list
            A list of sorted `int` chain ids of the given list of pdbs.
        """
        macromolecule_types = set(['Polyribonucleotide (RNA)','polyribonucleotide'])
        macromolecule_types.add('DNA/RNA Hybrid')
        macromolecule_types.add('NA-hybrid')
        macromolecule_types.add('polydeoxyribonucleotide/polyribonucleotide hybrid')
        macromolecule_types.add('Polydeoxyribonucleotide (DNA)')
        macromolecule_types.add('polydeoxyribonucleotide')

        with self.session() as session:
            if isinstance(pdbs, basestring):
                query = session.query(mod.ChainInfo.chain_name,
                     mod.ChainInfo.chain_id).\
                filter(mod.ChainInfo.entity_macromolecule_type.in_(macromolecule_types)).\
                filter_by(pdb_id=pdbs)
            else:
                query = session.query(mod.ChainInfo.chain_name,
                     mod.ChainInfo.chain_id).\
                filter(mod.ChainInfo.entity_macromolecule_type.in_(macromolecule_types)).\
                filter(mod.ChainInfo.pdb_id.in_(pdbs))

            chain_id = []
            for result in query:
                entry = result.chain_id
                chain_id.append(entry)
            return sorted(set(chain_id))

    def query(self, session, chain_id):
        """Create a query to find all mapped chains. This will produce a query
        to find all entries in ``exp_seq_chain_mappings``.

        Parameters
        ----------
        session : pymotifs.core.db.Session
            The session to use
        chain_id : int
            The chain id to use

        Returns
        -------
        query : Query
            The query.
        """
        return session.query(mod.ExpSeqChainMapping).\
            filter_by(chain_id=chain_id)

    def exp_id_old(self, chain_id):
        """Compute the experimetnal sequence id for the given chain id. This
        will look up all experimental sequences with the same sequence as the
        given chain id.

        Parameters
        ----------
        chain_id : int
            The chain id.

        Returns
        -------
        exp_seq_ids : list
            List of int experimental sequence ids
        """

        with self.session() as session:
            exp = mod.ExpSeqInfo
            query = session.query(exp.exp_seq_id).\
                join(mod.ChainInfo,
                     mod.ChainInfo.sequence == exp.sequence).\
                filter(mod.ChainInfo.chain_id == chain_id)

            if query.count() != 1:
                raise core.InvalidState("There should be exactly one matching"
                                        " experimental sequence")
            return query.one().exp_seq_id
        
    def exp_id_two_queries(self, chain_id):
        simplify_type = {}
        simplify_type['Polydeoxyribonucleotide (DNA)'] = 'dna'
        simplify_type['Polyribonucleotide (RNA)'] = 'rna'
        simplify_type['polyribonucleotide'] = 'rna'
        simplify_type['polydeoxyribonucleotide'] = 'dna'
        simplify_type['polydeoxyribonucleotide/polyribonucleotide hybrid'] = 'hybrid'
        simplify_type['DNA/RNA Hybrid'] = 'hybrid'
        data = []
        with self.session() as session:
            query = session.query(mod.ChainInfo.sequence,mod.ChainInfo.entity_macromolecule_type).\
                                    filter(mod.ChainInfo.chain_id == chain_id)
            for result in query:
                query2 = session.query(mod.ExpSeqInfo.exp_seq_id).\
                        filter(mod.ExpSeqInfo.sequence == result.sequence).filter(mod.ExpSeqInfo.entity_type == simplify_type[result.entity_macromolecule_type])
                if query2.count() != 1:
                    raise core.InvalidState("There should be exactly one matching"
                                            " experimental sequence")
                else:
                    data.append(result2.exp_seq_id for result2 in query2)
            return data  

    def exp_id(self, chain_id):
        """Compute the experimetnal sequence id for the given chain id. This
        will look up all experimental sequences with the same sequence and entity types
        as thegiven chain id.
        Parameters
        ----------
        chain_id : int
            The chain id.

        Returns
        -------
        exp_seq_ids : list
            List of int experimental sequence ids
        """
        simplify_type = {}
        simplify_type['Polydeoxyribonucleotide (DNA)'] = 'dna'
        simplify_type['Polyribonucleotide (RNA)'] = 'rna'
        simplify_type['polyribonucleotide'] = 'rna'
        simplify_type['polydeoxyribonucleotide'] = 'dna'
        simplify_type['polydeoxyribonucleotide/polyribonucleotide hybrid'] = 'hybrid'
        simplify_type['DNA/RNA Hybrid'] = 'hybrid'  

        with self.session() as session:
            exp = mod.ExpSeqInfo
            query = session.query(exp.exp_seq_id).\
                join(mod.ChainInfo,
                     # (mod.ChainInfo.sequence == exp.sequence) & (simplify_type[mod.ChainInfo.entity_macromolecule_type] == exp.entity_type)).\
                     (mod.ChainInfo.sequence == exp.sequence) & (case(simplify_type, value = mod.ChainInfo.entity_macromolecule_type) == exp.entity_type)).\
                filter(mod.ChainInfo.chain_id == chain_id)

            if query.count() != 1:
                raise core.InvalidState("There should be exactly one matching"
                                        " experimental sequence and entity type")
            return query.one().exp_seq_id            

    def data(self, chain_id, **kwargs):
        """Compute the mapping between the chain and an experimental sequence.

        Parameters
        ----------
        chain_id : int
            The chain id to use.

        Returns
        -------
        mapping : ExpSeqMapping
            A mapping object between the chain and the experimental sequence
            with the same sequence.
        """

        return mod.ExpSeqChainMapping(exp_seq_id=self.exp_id(chain_id),
                                      chain_id=chain_id)
