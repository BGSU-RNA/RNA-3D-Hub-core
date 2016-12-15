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
    pass


def as_mapped(result):
    return MappedChain(**row2dict(result))

DELETE = """
DELETE exp_seq_unit_mapping
FROM exp_seq_unit_mapping
JOIN unit_info
ON
    unit_info.unit_id = exp_seq_unit_mapping.unit_id
WHERE unit_info.pdb_id = :pdb_id
;
"""


class Loader(core.Loader):
    """The loader to use. This is a bit unusual in that we can't inherit from
    simple loader as we have to parse CIF files. Since this can take a long
    time for large files we don't process a single chain at a time, which would
    be easier to process. Instead we process a whole file at once.
    """

    dependencies = set([InfoLoader, PositionLoader, ChainMappingLoader,
                        UnitLoader])

    def has_data(self, pdb, **kwargs):
        """Checks if we have all the data for a given structure. This is a bit
        strange as compared to other has_data methods as it will raise skip for
        things without rna chains, such as unknown PDB ids. This could cause
        some issue with partially mapped structures, but I haven't seen this
        yet.

        Parameters
        ----------
        pdb : str
            The pdb id.

        Returns
        -------
        has : bool
            True if this has data for the given pdb.
        """

        with self.session() as session:
            query = session.query(mod.ChainInfo.chain_name).\
                filter(mod.ChainInfo.pdb_id == pdb).\
                filter(mod.ChainInfo.entity_macromolecule_type == 'Polyribonucleotide (RNA)').\
                distinct()

            possible = set(result.chain_name for result in query)

        if not possible:
            raise core.Skip("No chains to map")

        with self.session() as session:
            query = session.query(mod.UnitInfo.chain).\
                join(mod.ExpSeqUnitMapping,
                     mod.ExpSeqUnitMapping.unit_id == mod.UnitInfo.unit_id).\
                join(mod.ExpSeqPosition,
                     mod.ExpSeqPosition.exp_seq_position_id == mod.ExpSeqUnitMapping.exp_seq_position_id).\
                filter(mod.UnitInfo.pdb_id == pdb).\
                distinct()

            known = set(result.chain for result in query)

        return possible == known

    def remove(self, pdb, **kwargs):
        """Delete all mapping data for the given pdb.

        Parameters
        ----------
        pdb : str
            The pdb id to use.
        """

        with self.session() as session:
            session.execute(DELETE, {'pdb_id': pdb})

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

    def mapped_chains(self, pdb):
        """Get all mapped chains for the given pdb. This will look up all
        chain names that have been mapped to experimental sequences.

        Parameters
        ----------
        pdb : str
            The pdb id

        Returns
        -------
        chains : list
            A list of mapped chain names.
        """

        with self.session() as session:
            query = session.query(mod.ChainInfo.chain_name.label('name'),
                                  mod.ChainInfo.chain_id,
                                  mod.ExpSeqChainMapping.exp_seq_chain_mapping_id.label('id'),
                                  ).\
                join(mod.ExpSeqChainMapping,
                     mod.ExpSeqChainMapping.chain_id == mod.ChainInfo.chain_id).\
                filter(mod.ChainInfo.pdb_id == pdb).\
                filter(mod.ChainInfo.entity_macromolecule_type == 'Polyribonucleotide (RNA)')
            return sorted(as_mapped(result) for result in query)

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
