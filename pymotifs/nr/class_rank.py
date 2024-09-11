"""Store the assignments of chains to nr classes. This means the chains which
are part of an equivalence class. This requires that there is cached NR data to
store.
"""

import datetime as dt

from pymotifs import core

from pymotifs import models as mod
from pymotifs.constants import NR_CACHE_NAME

from pymotifs.nr.utils import BaseLoader
from pymotifs.nr.id_mapping import Loader as IdMappingLoader


class Loader(BaseLoader):

    dependencies = set([IdMappingLoader])
    # update_gap = dt.timedelta(7)  # Only update every 7 days
    allow_no_data = True
    mark = False

    @property
    def table(self):
        return mod.NrClassRank


    def has_data(self, *args, **kwargs):
        """
        We always want stages.py to think that data needs to be computed
        for this stage, so it does not get skipped.
        This works better than trying to manipulate the query method.
        """
        return False


    # def nr_class_name_checking(self):
    #     """
    #     Look for all existing nr_class_names in the nr_class_rank table.
    #     """
    #     with self.session() as session:
    #         query = session.query(mod.NrClassRank.nr_class_name).distinct()
    #     return [row.nr_class_name for row in query]


    def chains(self, grouping):
        """
        Compute the chain level data to store. The produced chains will be
        suitable for writing to the ``nr_class_rank`` table in the database.

        Parameters
        ----------
        release : str
            The nr release id.

        grouping : dict
            The dict with a 'members' entry that contains the chain to store.

        Returns
        -------
        chains : list
            A list of chain dicts to store.
        """

        data = []
        for group in grouping:
            for chain in group['members']:
                data.append({
                    'ife_id': chain['id'],
                    'nr_class_name': group['name']['full'],
                    'rank': chain['rank'],
                    'nr_class_id': group['name']['class_id']
                })
        return data


    def data(self, release, **kwargs):
        """
        Compute the data to store.

        Parameters
        ----------
        release : str
            The nr release id to use.

        Returns
        -------
        chains : list
            A list of chain dicts to store.
        """

        data = self.cached(NR_CACHE_NAME)
        if not data:
            raise core.Skip("Nothing to do here, maybe too few files, maybe an earlier stage failed")
        return self.chains(data['groups'])
