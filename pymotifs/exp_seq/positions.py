from pymotifs import core

from pymotifs.models import ExpSeqInfo as Exp
from pymotifs.models import ExpSeqPosition as Position
from pymotifs.models import ExpSeqChainMapping as Mapping
from pymotifs.models import ChainInfo

from pymotifs.exp_seq.info import Loader as InfoLoader
from pymotifs.exp_seq.chain_mapping import Loader as ExpMappingLoader


class Loader(core.Loader):
    dependencies = set([ExpMappingLoader, InfoLoader])

    def to_process(self, pdbs, **kwargs):
        with self.session() as session:
            query = session.query(Mapping.exp_seq_id).\
                join(ChainInfo, ChainInfo.id == Mapping.chain_id).\
                filter(ChainInfo.pdb_id.in_(pdbs)).\
                distinct()

            return [result.exp_seq_id for result in query]

    # def known(self, pdb):
    #     with self.session() as session:
    #         query = session.query(Position.exp_seq_id).\
    #             join(Mapping, Position.exp_seq_id == Mapping.exp_seq_id).\
    #             join(ChainInfo, ChainInfo.id == Mapping.chain_id).\
    #             filter(ChainInfo.pdb_id == pdb).\
    #             distinct()

    #     return set(result.exp_seq_id for result in query)

#     def missing(self, pdb):
#         helper = Structure(self.session.maker)
#         possible = set(p[1] for p in helper.rna_chains(pdb, return_id=True))
#         return possible - set(self.known(pdb))

    def has_data(self, exp_seq_id, **kwargs):
        with self.session() as session:
            query = session.query(Position).\
                filter(Position.exp_seq_id == exp_seq_id)

            return bool(query.count())

    def remove(self, exp_seq_id, **kwargs):
        with self.session() as session:
            session.query(Position).\
                filter(Position.exp_seq_id == exp_seq_id).\
                delete(synchronize_session=False)

    def sequence(self, exp_seq_id):
        with self.session() as session:
            return session.query(Exp).get(exp_seq_id).sequence

    def positions(self, exp_id, sequence):
        positions = []
        for index, char in enumerate(sequence):
            positions.append({
                'exp_seq_id': exp_id,
                'unit': char,
                'index': index
            })
        return positions

    def data(self, exp_seq_id, **kwargs):
        data = []
        sequence = self.sequence(exp_seq_id)
        for position in self.positions(exp_seq_id, sequence):
            data.append(Position(**position))

        return data
