"""Load data about each experimental sequence positions. This will process all
experimental sequences and write data on each position in the sequence.
"""

from pymotifs import core
from pymotifs import models as mod

from pymotifs.exp_seq.info import Loader as InfoLoader
from pymotifs.exp_seq.chain_mapping import Loader as ExpMappingLoader


class Loader(core.SimpleLoader):
    dependencies = set([ExpMappingLoader, InfoLoader])

    @property
    def table(self):
        return mod.ExpSeqPosition

    mark = False

    def to_process(self, pdbs, **kwargs):
        """Find all stored experimental sequences for the given pdb ids. This
        will only use experimental sequences which show up in the given
        structures.

        Parameters
        ----------
        pdbs : list
            List of pdb ids to use

        Returns
        -------
        exp_seq_ids : list
            A list of experimental sequence ids to process.
        """

        with self.session() as session:
            query = session.query(mod.ExpSeqChainMapping.exp_seq_id).\
                join(mod.ChainInfo, mod.ChainInfo.chain_id == mod.ExpSeqChainMapping.chain_id).\
                filter(mod.ChainInfo.pdb_id.in_(pdbs)).\
                distinct()

            return [result.exp_seq_id for result in query]

    def query(self, session, exp_seq_id):
        """Compute a query for all store experimental positions with the given
        experimetnal sequence id.

        Parameters
        ----------
        session : pymotifs.core.db.Session
            The session wrapper to use.
        exp_seq_id : int
            The experimental sequence id

        Returns
        -------
        query : Query
            A query for all experimental sequence positions for the given id.
        """
        return session.query(mod.ExpSeqPosition).\
            filter(mod.ExpSeqPosition.exp_seq_id == exp_seq_id)

    def sequence(self, exp_seq_id):
        """Get the sequence for the given experimental sequence id. This will
        raise an error if the exp_seq_id does not exist in the database.

        Parameters
        ----------
        exp_seq_id : int
            The id.

        Returns
        -------
        seq : str
            The sequence
        """
        with self.session() as session:
            exp = session.query(mod.ExpSeqInfo).get(exp_seq_id)
            return exp.sequence

    def positions(self, exp_id, sequence):
        """Compute a dictonary for each position in an experimental sequence.
        Each dictionary in the resulting list will have the a 'exp_seq_id',
        'unit' (the character at that position', 'normalized_unit' (the
        normalized unit), 'index' (index of the position). The normalization is
        done using the normalization from
        `pymotifs.exp_seq.info.Loader.translate`.

        Parameters
        ----------
        exp_id : int
            The experimental sequence id.
        sequence : str
            The sequence of the experimental sequence.

        Returns
        -------
        positions : list
            A list of dicts.
        """

        positions = []
        info = InfoLoader(self.config, self.session)
        for index, char in enumerate(sequence):
            norm_char = info.translate(char)

            positions.append({
                'exp_seq_id': exp_id,
                'unit': char,
                'normalized_unit': norm_char,
                'index': index
            })
        return positions

    def data(self, exp_seq_id, **kwargs):
        """Compute the positions for the given experimental sequences. This
        returns the same data as `Loader.positions`.

        Returns
        -------
        positions : list
            A list of dictionaries as from `Loader.positions`.
        """

        sequence = self.sequence(exp_seq_id)
        return self.positions(exp_seq_id, sequence)
