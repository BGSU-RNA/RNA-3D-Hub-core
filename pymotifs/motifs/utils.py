"""This is a module that contains the basic logic for import motif data to the
pipeline. It does not contain any runnable stages, but it is used for the other
motif importers.
"""

import abc

from pymotifs import core
from pymotifs import models as mod

from pymotifs.motifs.release import Loader as ReleaseLoader


class BaseLoader(core.SimpleLoader):
    """This is a convience class for the motif atlas loaders other than the
    release loader and cleanup. They all have the same logic to them for
    detecting if we need to update and what to process. Thus we create a base
    class which handles that logic.

    All instances of this class MUST have a a table property which is the table
    to write to as this is used in the query method.
    """

    __metaclass__ = abc.ABCMeta

    def mark_processed(self, pair, **kwargs):
        super(BaseLoader, self).mark_processed(pair[0], **kwargs)

    def to_process(self, pdbs, **kwargs):
        """Compute the data to process. The input PDB's are ignored and instead
        the cache is examine for motif data, IL and HL to import. This will
        produce a list of tuples of the motifs to import.

        Parameters
        ----------
        pdbs : list
            Ingored

        Returns
        -------
        A list of tuples like [('IL', '1.0'), ('HL', '1.0'] to import.
        """
        current, _ = ReleaseLoader(self.config, self.session).current_id()
        data = []
        for loop_type in ReleaseLoader.types:
            cached = self.cached(loop_type)
            if not cached:
                raise core.InvalidState("No cached data")

            if cached['release'] != current:
                raise core.InvalidState("Caching does not match excepted ID")
            data.append((loop_type, current))
        return data

    def __self_query__(self, session, pair):
        """Create a query against the table itself to import.
        """
        info = mod.MlMotifsInfo
        return session.query(info).\
            filter(info.type == pair[0]).\
            filter(info.ml_release_id == pair[1])

    def query(self, session, pair):
        """Create a query to check if this has already import data. This uses
        the loop_type and the release to check if data has already been
        processed. It assumes there are motif_id and ml_release_id
        columns which are used to check if we have entries of the right data.

        Parameters
        ----------
        session : pymotifs.core.db.Session
            The session to use

        pair : tuple
            A tuple produced by to_process

        Returns
        -------
        A query that checks if there is already data for this tuple in the
        given table.
        """
        info = mod.MlMotifsInfo
        if self.table == info:
            return self.__self_query__(session, pair)

        loop_type, release = pair
        return session.query(self.table).\
            filter(self.table.motif_id.like(loop_type + '_%')).\
            filter(self.table.ml_release_id == release)

    @abc.abstractmethod
    def data(self, pair, **kwargs):
        pass
