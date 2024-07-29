"""
Map experimental sequence positions to unit ids.
This will process the pdbx_poly_seq_scheme entry in cif files
to produce a mapping between experimental sequence positions and unit ids.
It deals with positions that are not mapped to unit ids as well.
"""

from collections import defaultdict
from collections import namedtuple as nt
from sqlalchemy import and_

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import row2dict

from pymotifs.exp_seq.info import Loader as InfoLoader
from pymotifs.exp_seq.positions import Loader as PositionLoader
from pymotifs.exp_seq.chain_mapping import Loader as ChainMappingLoader
from pymotifs.units.info import Loader as UnitLoader

nucleic_acid_types = set(['Polyribonucleotide (RNA)','polyribonucleotide'])
nucleic_acid_types.add('DNA/RNA Hybrid')
nucleic_acid_types.add('NA-hybrid')
nucleic_acid_types.add('polydeoxyribonucleotide/polyribonucleotide hybrid')
nucleic_acid_types.add('Polydeoxyribonucleotide (DNA)')
nucleic_acid_types.add('polydeoxyribonucleotide')


class MappedChain(nt('MappedChain', ['id', 'chain_id', 'name'])):
    """
    id is exp_seq_chain_mapping_id
    chain_id is chain_id
    name is chain_name
    """
    @classmethod
    def from_dict(cls, result):
        return cls(**row2dict(result))
    pass


class Loader(core.SimpleLoader):
    """
    The loader to use.
    This is a bit unusual in that we can't inherit from
    simple loader as we have to parse CIF files.
    Since this can take a long time for large files we don't process
    a single chain at a time, which would be easier to process.
    Instead we process a whole file at once.

    However, the to_process method identifies individual chains
    that need to be mapped.
    """

    dependencies = set([InfoLoader, PositionLoader, ChainMappingLoader,
                        UnitLoader])

    mark = False  # OK for stages where to_process returns a list of pdb ids,
                  # but this stage works at the chain level; can't use this


    def to_process(self, pdbs, **kwargs):
        """
        Make a list of tuples of pdb ids and chains that are not already mapped
        This list will be processed one by one by the data method.
        Example: ('4V9F',['0','9'])
        """

        # find all rna/dna chains in the given pdbs that are already in exp_seq_pdb
        # query takes less than 1 second
        with self.session() as session:
            query = session.query(mod.ExpSeqPdb.pdb_id, mod.ExpSeqPdb.chain_name).\
                filter(mod.ExpSeqPdb.pdb_id.in_(pdbs))
            pdb_id_to_chains = defaultdict(set)
            for result in query:
                pdb_id_to_chains[result.pdb_id].add(result.chain_name)

        self.logger.info('Found %d pdb ids with chains in exp_seq_pdb' % len(pdb_id_to_chains))

        # find all chains in the given pdbs that already have unit mappings
        # query takes about 66 seconds for 8000 pdbs
        with self.session() as session:
            query = session.query(mod.ExpSeqPdb.pdb_id, mod.ExpSeqPdb.chain_name).\
                join(mod.ExpSeqChainMapping, mod.ExpSeqPdb.chain_id == mod.ExpSeqChainMapping.chain_id).\
                join(mod.ExpSeqUnitMapping, mod.ExpSeqChainMapping.exp_seq_chain_mapping_id == mod.ExpSeqUnitMapping.exp_seq_chain_mapping_id).\
                filter(mod.ExpSeqPdb.pdb_id.in_(pdbs))
            pdb_id_to_mapped_chains = defaultdict(set)
            for result in query:
                pdb_id_to_mapped_chains[result.pdb_id].add(result.chain_name)

        self.logger.info('Found %d pdb ids with mapped chains' % len(pdb_id_to_mapped_chains))

        # identify the pdbs that have chains that are not mapped
        pdb_id_to_tentative_chains = {}
        tentative_pdbs = []
        for pdb_id in pdbs:
            # self.logger.info('PDB %s has exp_seq chains %s' % (pdb, sorted(pdb_id_to_chains[pdb])))
            # self.logger.info('PDB %s has exp_seq_unit_mapping chains %s' % (pdb, sorted(pdb_id_to_mapped_chains[pdb])))
            chains_to_process = pdb_id_to_chains[pdb_id] - pdb_id_to_mapped_chains[pdb_id]
            if chains_to_process:
                pdb_id_to_tentative_chains[pdb_id] = chains_to_process
                tentative_pdbs.append(pdb_id)
                self.logger.info('PDB %s has chains %s to process' % (pdb_id, sorted(chains_to_process)))

        if len(tentative_pdbs) == 0:
            raise core.Skip("All chains already mapped")

        # double check that the pdbs we have, actually have entries in unit_info
        # Strange cases like 7PJS chains v, w, z are not in unit_info
        # this query returns chains we could actually map
        # this query is slow so we only do it on the tentative pdbs
        with self.session() as session:
            query = session.query(mod.ExpSeqPdb.chain_name, mod.ExpSeqPdb.pdb_id).\
                join(mod.UnitInfo, and_(mod.UnitInfo.pdb_id == mod.ExpSeqPdb.pdb_id, mod.UnitInfo.chain == mod.ExpSeqPdb.chain_name)).\
                filter(mod.ExpSeqPdb.pdb_id.in_(tentative_pdbs))
            pdb_id_to_confirmed_chains = defaultdict(set)
            for result in query:
                pdb_id_to_confirmed_chains[result.pdb_id].add(result.chain_name)

        all_pdb_chains = []
        for pdb_id in tentative_pdbs:
            # self.logger.info('PDB %s has exp_seq chains %s' % (pdb, sorted(pdb_id_to_chains[pdb])))
            # self.logger.info('PDB %s has exp_seq_unit_mapping chains %s' % (pdb, sorted(pdb_id_to_mapped_chains[pdb])))
            chains_to_process = pdb_id_to_tentative_chains[pdb_id] & pdb_id_to_confirmed_chains[pdb_id]
            if chains_to_process:
                all_pdb_chains.append((pdb_id, sorted(chains_to_process)))
                self.logger.info('PDB %s has chains %s that can be processed' % (pdb_id, sorted(chains_to_process)))
            else:
                self.logger.info('PDB %s has chains %s that cannot be processed' % (pdb_id, sorted(pdb_id_to_tentative_chains[pdb_id])))

        if len(all_pdb_chains) == 0:
            raise core.Skip("All chains already mapped")

        return all_pdb_chains


    def query(self, session, pdb_chains):
        """
        Query the database for mappings for the given (PDB,chain_names) tuple.
        This way, the new to_process method can identify individual chains
        which did not get mapped.
        For some reason, that happens, and then those chains never get into
        the exp_seq_unit_mapping table.

        Parameters
        ----------
        session : pymotifs.core.db.Session
            The session to use.
        pdb_chains : tuple of (pdb_id, chain_names)
            pdb_id : string, The PDB id to use
            chain_names : list of strings, The chain names to use

        Returns
        -------
        A query that will find mappings for the given structure.
        """

        pdb, chain_names = pdb_chains

        return session.query(mod.ExpSeqUnitMapping).\
            join(mod.ExpSeqChainMapping,
                 mod.ExpSeqChainMapping.exp_seq_chain_mapping_id == mod.ExpSeqUnitMapping.exp_seq_chain_mapping_id).\
            join(mod.ChainInfo,
                 mod.ChainInfo.chain_id == mod.ExpSeqChainMapping.chain_id).\
            filter(mod.ChainInfo.pdb_id == pdb).\
            filter(mod.ChainInfo.chain_name.in_(chain_names))


    def mapped_chains(self, pdb, chain_names, extended=True):
        """
        Get all desired mapped chains for the given pdb. This will look up all
        chain names that have been mapped to experimental sequences.

        Parameters
        ----------
        pdb : str
            The pdb id
        extended : bool
            Work with only "Polyribonucleotide (RNA)" chains (False) or
            allow additional chains to be considered (True; default).

        Returns
        -------
        chains : list
            A list of lists containing 'id', 'chain_id', 'name'
        """

        with self.session() as session:
            query = session.query(mod.ChainInfo.chain_name.label('name'),
                                  mod.ChainInfo.chain_id,
                                  mod.ExpSeqChainMapping.exp_seq_chain_mapping_id.label('id'),
                                  ).\
                join(mod.ExpSeqChainMapping,
                     mod.ExpSeqChainMapping.chain_id == mod.ChainInfo.chain_id).\
                filter(mod.ChainInfo.pdb_id == pdb).\
                filter(mod.ChainInfo.chain_name.in_(chain_names)).\
                filter(mod.ChainInfo.entity_macromolecule_type.in_(nucleic_acid_types))
            return sorted(MappedChain.from_dict(result) for result in query)


    def exp_mapping(self, pdb, chains):
        """
        Compute a mapping from index in a chain to experimental sequence
        position.
        This will look up all positions and ...

        Parameters
        ----------
        pdb : str
            Pdb id to use
        chains : list
            List of chain names.

        Returns
        -------
        mapping : dict
            A dictonary of mappings from (chain_name, index) to experimental
            sequence id.
        """

        chain_names = chains
        with self.session() as session:
            query = session.query(mod.ExpSeqPosition.index,
                                  mod.ExpSeqPosition.exp_seq_position_id,
                                  mod.ChainInfo.chain_name).\
                join(mod.ExpSeqChainMapping,
                     mod.ExpSeqChainMapping.exp_seq_id == mod.ExpSeqPosition.exp_seq_id).\
                join(mod.ChainInfo,
                     mod.ChainInfo.chain_id == mod.ExpSeqChainMapping.chain_id).\
                filter(mod.ChainInfo.pdb_id == pdb).\
                filter(mod.ChainInfo.chain_name.in_(chain_names))

            seen = set()
            chain_index_to_exp_seq_position_id = {}
            for result in query:
                # key is chain_name like 'A' and index runs from 0 to length-1
                # add 1 to this index so that it is the same as chain index in unit_info
                key = (result.chain_name, result.index + 1)
                chain_index_to_exp_seq_position_id[key] = result.exp_seq_position_id
                seen.add(result.chain_name)

            if seen != set(chain_names):
                msg = "Could not get mappings for all chains in %s, %s"
                raise core.InvalidState(msg, pdb, ', '.join(chain_names))

            return chain_index_to_exp_seq_position_id


    def chain_mapping(self, pdb, mapped_chains, exp_mapping):
        """
        Compute the mapping between experimental sequence position id and
        unit id.

        2024-06-19 CLZ: We should be able to get the same index numbers from
        the unit_info table.  They got there from the .cif file in the first place.

        Parameters
        ----------
        cif : fr3d.cif.reader.Cif
            The cif data structure to use.
        chain : list
            List of chain names to look up.
        exp_seq_mapping :
            A dict from `exp_mapping` that maps to experimental sequence
            position id.

        Yields
        -------
        mapping : ExpSeqUnitMapping
            A series of unit mappings that map from experimental sequence
            position to unit id.
        """

        chains = []
        for mapped_chain in mapped_chains:
            chains.append(mapped_chain.name)

        # new code to get the chain_index from unit_info and avoid reading the cif file
        with self.session() as session:
            query = session.query(mod.UnitInfo.unit_id,mod.UnitInfo.chain,mod.UnitInfo.chain_index).\
                filter(mod.UnitInfo.pdb_id == pdb).\
                filter(mod.UnitInfo.chain.in_(chains)).\
                order_by(mod.UnitInfo.chain,mod.UnitInfo.chain_index)

            chain_to_index_to_unit_id = {}
            for result in query:
                if result.chain_index is not None:
                    if not result.chain in chain_to_index_to_unit_id:
                        chain_to_index_to_unit_id[result.chain] = {}
                    chain_to_index_to_unit_id[result.chain][result.chain_index] = result.unit_id

                # self.logger.info('unit_info: %s %s %s' % (result.chain, result.unit_id, result.chain_index))

        for mapped_chain in mapped_chains:
            chain = mapped_chain.name
            exp_seq_chain_mapping_id = mapped_chain.id
            for index, unit_id in chain_to_index_to_unit_id[chain].items():

                if (chain, index) not in exp_mapping:
                    raise core.InvalidState("No exp_seq_position_id for chain %s index %s" % (chain, index))

                # next line is the hard part; this is an integer that is unique to the
                # chain and the position in the chain, over all NA chains in the database
                exp_seq_position_id = exp_mapping[(chain, index)]

                self.logger.info('Mapping %s to %s in chain %s with exp_seq_chain_mapping_id %s' % (unit_id, exp_seq_position_id, chain, exp_seq_chain_mapping_id))

                yield mod.ExpSeqUnitMapping(
                    unit_id=unit_id,
                    exp_seq_chain_mapping_id=exp_seq_chain_mapping_id,
                    exp_seq_position_id=exp_seq_position_id,
                    chain=chain,
                    )


    def data(self, pdb_chains, **kwargs):
        """
        Compute the data for the given tuple of pdb id and list of chain names.
        This will load the cif file and
        get the mapping for residues in all designated NA chains.

        Yields
        ------
        entry : ExpSeqUnitMapping
            An entry as from `Loader.chain_mapping`.
        """

        pdb, chain_names = pdb_chains

        # old code below checked to see if the chains were registered, but did
        # not check to see if the chains were already mapped.
        # Now we know which chains actually need to be mapped in this stage.

        # look up chain id and id
        chains = self.mapped_chains(pdb, chain_names)

        if not chains:
            self.logger.info('No chains found in %s' % pdb)
            raise core.InvalidState("No chains found in %s" % pdb)

        self.logger.info('chain tuples from mapped_chains: %s' % str(chains))

        exp_mapping = self.exp_mapping(pdb, chain_names)

        self.logger.info('exp_seq/mapping exp_mapping:')
        self.logger.info(sorted(exp_mapping.items()))

        for entry in self.chain_mapping(pdb, chains, exp_mapping):
            yield entry
