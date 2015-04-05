from pymotifs import core

from pymotifs.models import ChainInfo
from pymotifs.models import UnitInfo
from pymotifs.models import ExpSeqPosition as ExpPosition
from pymotifs.models import ExpSeqUnitMapping as UnitMapping
from pymotifs.models import ExpSeqChainMapping as ChainMapping

from pymotifs.exp_seq.info import Loader as InfoLoader
from pymotifs.units.info import Loader as UnitLoader


class Loader(core.Loader):
    insert_max = 5000
    dependencies = set([InfoLoader, UnitLoader])

    def has_data(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(UnitMapping).\
                join(UnitInfo, UnitInfo.id == UnitMapping.unit_id).\
                filter(UnitInfo.pdb_id == pdb)
            return bool(query.count())

    def remove(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(UnitInfo.id).filter(UnitInfo.pdb_id == pdb)
            units = [result.id for result in query]

        with self.session() as session:
            session.query(UnitMapping).\
                filter(UnitMapping.unit_id.in_(units)).\
                delete(synchronize_session=False)

    def chain_mapping(self, cif, chain, exp_mapping):
        units = []
        seq = cif.experimental_sequence_mapping(chain)
        for index, (_, _, unit_id) in enumerate(seq):
            if not unit_id:
                unit_id = None
            seq_id = exp_mapping[index]
            units.append(UnitMapping(unit_id=unit_id,
                                     exp_seq_position_id=seq_id))
        return units

    def exp_mapping(self, pdb, chain):
        with self.session() as session:
            query = session.query(ExpPosition.index, ExpPosition.id).\
                join(ChainMapping,
                     ChainMapping.exp_seq_id == ExpPosition.exp_seq_id).\
                join(ChainInfo, ChainInfo.id == ChainMapping.chain_id).\
                filter(ChainInfo.pdb_id == pdb).\
                filter(ChainInfo.chain_name == chain)

            mapping = {}
            for result in query:
                mapping[result.index] = result.id
            return mapping

    def mapped_chains(self, pdb):
        with self.session() as session:
            query = session.query(ChainInfo.chain_name).\
                join(ChainMapping, ChainMapping.chain_id == ChainInfo.id).\
                filter(ChainInfo.pdb_id == pdb)
            return [result.chain_name for result in query]

    def data(self, pdb, **kwargs):
        cif = self.cif(pdb)
        chains = self.mapped_chains(cif.pdb)

        data = []
        for chain_name in chains:
            exp_mapping = self.exp_mapping(cif.pdb, chain_name)
            mapping = self.chain_mapping(cif, chain_name, exp_mapping)
            data.extend(mapping)
        return data
