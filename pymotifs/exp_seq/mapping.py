"""Map experimental sequence positions to unit ids. This will process the
pdbx_poly_seq_scheme entry in cif files to produce a mapping between
experimental sequence positions and unit ids. It deals with positions that are
not mapped to unit ids as well.
"""

from collections import namedtuple as nt

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import row2dict

from pymotifs.exp_seq.info import Loader as InfoLoader
from pymotifs.exp_seq.positions import Loader as PositionLoader
from pymotifs.exp_seq.chain_mapping import Loader as ChainMappingLoader
from pymotifs.units.info import Loader as UnitLoader


class MappedChain(nt('MappedChain', ['id', 'chain_id', 'name'])):
    @classmethod
    def from_dict(cls, result):
        return cls(**row2dict(result))
    pass


class Loader(core.SimpleLoader):
    """The loader to use. This is a bit unusual in that we can't inherit from
    simple loader as we have to parse CIF files. Since this can take a long
    time for large files we don't process a single chain at a time, which would
    be easier to process. Instead we process a whole file at once.
    """

    dependencies = set([InfoLoader, PositionLoader, ChainMappingLoader,
                        UnitLoader])

    def query(self, session, pdb):
        """Query the database for mappings for the given PDB.

        Parameters
        ----------
        session : pymotifs.core.db.Session
            The session to use.
        pdb : str
            The PDB id to use

        Returns
        -------
        A query that will find mappings for the given structure.
        """
        return session.query(mod.ExpSeqUnitMapping).\
            join(mod.ExpSeqChainMapping,
                 mod.ExpSeqChainMapping.exp_seq_chain_mapping_id == mod.ExpSeqUnitMapping.exp_seq_chain_mapping_id).\
            join(mod.ChainInfo,
                 mod.ChainInfo.chain_id == mod.ExpSeqChainMapping.chain_id).\
            filter(mod.ChainInfo.pdb_id == pdb)

    def chain_mapping(self, cif, mapped_chains, exp_mapping):
        """Compute the mapping between experimental sequence position id and
        unit id.

        Parameters
        ----------
        cif : fr3d.cif.reader.Cif
            The cif data structure to use.
        chain : list
            List of chain names to look up.
        exp_seq_mapping :
            A dict from `exp_mapping` that maps to experimental sequence
            position id.

        Yeilds
        -------
        mapping : ExpSeqUnitMapping
            A series of unit mappings that map from experimental sequence
            position to unit id.
        """

        trans = {m.name: m for m in mapped_chains}
        for mapping in cif.experimental_sequence_mapping(trans.keys()):
            unit_id = mapping['unit_id']
            index = mapping['index']
            chain = mapping['chain']
            key = (chain, index)

            if key not in exp_mapping:
                raise core.InvalidState("No pos id for %s" % str(key))

            pos_id = exp_mapping[key]
            mapped = trans[chain]
            yield mod.ExpSeqUnitMapping(
                unit_id=unit_id,
                exp_seq_chain_mapping_id=mapped.id,
                exp_seq_position_id=pos_id,
                chain=chain,
            )

    def exp_mapping(self, pdb, chains):
        """Compute a mapping from index in a chain to experimental sequence
        position. This will look up all positions and

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

        chain_names = [c.name for c in chains]
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
            mapping = {}
            for result in query:
                key = (result.chain_name, result.index)
                mapping[key] = result.exp_seq_position_id
                seen.add(result.chain_name)

            if seen != set(chain_names):
                msg = "Could not get mappings for all chains in %s, %s"
                raise core.InvalidState(msg, pdb, ', '.join(chain_names))

            return mapping

    def mapped_chains(self, pdb, extended=True):
        """Get all mapped chains for the given pdb. This will look up all
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
            A list of mapped chain names.
        """

        macromolecule_types = set(['Polyribonucleotide (RNA)'])
        if extended:
            macromolecule_types.add('DNA/RNA Hybrid')

        with self.session() as session:
            query = session.query(mod.ChainInfo.chain_name.label('name'),
                                  mod.ChainInfo.chain_id,
                                  mod.ExpSeqChainMapping.exp_seq_chain_mapping_id.label('id'),
                                  ).\
                join(mod.ExpSeqChainMapping,
                     mod.ExpSeqChainMapping.chain_id == mod.ChainInfo.chain_id).\
                filter(mod.ChainInfo.pdb_id == pdb).\
                filter(mod.ChainInfo.entity_macromolecule_type.in_(macromolecule_types))
            return sorted(MappedChain.from_dict(result) for result in query)

    def data(self, pdb, **kwargs):
        """Compute the data for the given pdb. This will load the cif file and
        get the mapping for residues in all RNA chains.

        Yields
        ------
        entry : ExpSeqUnitMapping
            An entry as from `Loader.chain_mapping`.
        """

        cif = self.cif(pdb)
        chains = self.mapped_chains(cif.pdb)
        if not chains:
            raise core.InvalidState("Found no chains in %s", pdb)
        exp_mapping = self.exp_mapping(cif.pdb, chains)
        for entry in self.chain_mapping(cif, chains, exp_mapping):
            yield entry
