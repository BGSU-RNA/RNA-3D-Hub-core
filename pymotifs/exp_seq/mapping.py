from pymotifs import core

from pymotifs.models import ChainInfo
from pymotifs.models import UnitInfo
from pymotifs.models import ExpSeqPosition as ExpPosition
from pymotifs.models import ExpSeqUnitMapping as UnitMapping
from pymotifs.models import ExpSeqChainMapping as ChainMapping

from pymotifs.exp_seq.info import Loader as InfoLoader
from pymotifs.exp_seq.positions import Loader as PositionLoader
from pymotifs.units.info import Loader as UnitLoader


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
            query = session.query(UnitInfo.unit_id).filter(UnitInfo.pdb_id == pdb)
            units = [result.id for result in query]

        with self.session() as session:
            session.query(UnitMapping).\
                filter(UnitMapping.unit_id.in_(units)).\
                delete(synchronize_session=False)

    def chain_mapping(self, cif, chain, exp_mapping):
        seq = cif.experimental_sequence_mapping(chain)

        for index, (_, _, unit_id) in enumerate(seq):
            if not unit_id:
                unit_id = None

            if long(index) not in exp_mapping:
                raise core.InvalidState("Could not find seq id for index %s" %
                                        index)

            seq_id = exp_mapping[long(index)]
            yield UnitMapping(unit_id=unit_id, exp_seq_position_id=seq_id)

    def exp_mapping(self, pdb, chain):
        with self.session() as session:
            query = session.query(ExpPosition.index, ExpPosition.exp_seq_position_id).\
                join(ChainMapping,
                     ChainMapping.exp_seq_id == ExpPosition.exp_seq_id).\
                join(ChainInfo, ChainInfo.chain_id == ChainMapping.chain_id).\
                filter(ChainInfo.pdb_id == pdb).\
                filter(ChainInfo.chain_name == chain)

            mapping = {}
            for result in query:
                mapping[result.index] = result.id
            return mapping

    def mapped_chains(self, pdb):
        with self.session() as session:
            query = session.query(ChainInfo.chain_name).\
                join(ChainMapping, ChainMapping.chain_id == ChainInfo.chain_id).\
                filter(ChainInfo.pdb_id == pdb).\
                filter(ChainInfo.entity_macromolecule_type == 'Polyribonucleotide (RNA)')
            return [result.chain_name for result in query]

    def data(self, pdb, **kwargs):

        cif = self.cif(pdb)
        chains = self.mapped_chains(cif.pdb)

        for chain_name in chains:
            exp_mapping = self.exp_mapping(cif.pdb, chain_name)
            if not exp_mapping:
                raise core.InvalidState("No mapping generated for %s, %s "
                                        "indicting no positions, skipping" %
                                        (cif.pdb, chain_name))

            for entry in self.chain_mapping(cif, chain_name, exp_mapping):
                yield entry
