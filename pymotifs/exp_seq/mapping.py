from pymotifs import core

from pymotifs.models import ChainInfo
from pymotifs.models import UnitInfo
from pymotifs.models import ExpSeqPosition as ExpPosition
from pymotifs.models import ExpSeqUnitMapping as UnitMapping
from pymotifs.models import ExpSeqChainMapping as ChainMapping

from pymotifs.exp_seq.info import Loader as InfoLoader
from pymotifs.exp_seq.positions import Loader as PositionLoader
from pymotifs.units.info import Loader as UnitLoader


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
    dependencies = set([InfoLoader, PositionLoader, UnitLoader])

    def has_data(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(ChainInfo.chain_name).\
                filter(ChainInfo.pdb_id == pdb).\
                filter(ChainInfo.entity_macromolecule_type == 'Polyribonucleotide (RNA)').\
                distinct()

            possible = set(result.chain_name for result in query)

        if not possible:
            raise core.Skip("No chains to map")

        with self.session() as session:
            query = session.query(UnitInfo.chain).\
                join(UnitMapping, UnitMapping.unit_id == UnitInfo.unit_id).\
                join(ExpPosition,
                     ExpPosition.exp_seq_position_id == UnitMapping.exp_seq_position_id).\
                filter(UnitInfo.pdb_id == pdb).\
                distinct()

            known = set(result.chain for result in query)

        return possible == known

    def remove(self, pdb, **kwargs):
        with self.session() as session:
            session.execute(DELETE, {'pdb_id': pdb})

    def chain_mapping(self, cif, chain, exp_mapping):
        for (_, seq_id, unit_id) in cif.experimental_sequence_mapping(chain):
            if not unit_id:
                unit_id = None

            parts = seq_id.split("|")
            index = long(parts[4]) - 1

            if index not in exp_mapping:
                raise core.InvalidState("No pos id for %s" % seq_id)

            pos_id = exp_mapping[index]
            yield UnitMapping(unit_id=unit_id, exp_seq_position_id=pos_id)

    def exp_mapping(self, pdb, chains):
        with self.session() as session:
            query = session.query(ExpPosition.index,
                                  ExpPosition.exp_seq_position_id,
                                  ChainInfo.chain_name).\
                join(ChainMapping,
                     ChainMapping.exp_seq_id == ExpPosition.exp_seq_id).\
                join(ChainInfo, ChainInfo.chain_id == ChainMapping.chain_id).\
                filter(ChainInfo.pdb_id == pdb).\
                filter(ChainInfo.chain_name.in_(chains))

            seen = set()
            mapping = {}
            for result in query:
                mapping[result.index] = result.exp_seq_position_id
                seen.add(result.chain_name)

            if seen != set(chains):
                msg = "Could not get mappings for all chains in %s, %s"
                raise core.InvalidState(msg, pdb, ', '.join(chains))

            return mapping

    def mapped_chains(self, pdb):
        with self.session() as session:
            query = session.query(ChainInfo.chain_name).\
                join(ChainMapping,
                     ChainMapping.chain_id == ChainInfo.chain_id).\
                filter(ChainInfo.pdb_id == pdb).\
                filter(ChainInfo.entity_macromolecule_type == 'Polyribonucleotide (RNA)')
            return [result.chain_name for result in query]

    def data(self, pdb, **kwargs):
        cif = self.cif(pdb)
        chains = self.mapped_chains(cif.pdb)
        if not chains:
            raise core.InvalidState("Found no chains in %s", pdb)
        exp_mapping = self.exp_mapping(cif.pdb, chains)
        print(exp_mapping)
        for entry in self.chain_mapping(cif, chains, exp_mapping):
            print(entry)
            yield entry
