"""Load the motif parent data. This will load from the cached data and then use
that create the parent level mappings. If there are no parents, which should be
rare, or if this first release this will skip processing.
"""

from pymotifs import core
from pymotifs import models as mod

from pymotifs.motifs.utils import BaseLoader
from pymotifs.motifs.info import Loader as InfoLoader
from pymotifs.motifs.release import Loader as ReleaseLoader


class Loader(BaseLoader):
    """Loader for the motif parent data.
    """
    dependencies = set([ReleaseLoader, InfoLoader])

    @property
    def table(self):
        return mod.MlParents

    def parents(self, cached):
        """Compute the parents for all cached motifs.

        Parameters
        ----------
        cached : dict
            The dictonary of cached data

        Returns
        --------
        parents : list
            List of parents dict. These dicts are suitable for storing into the
            ml_parents table.
        """

        data = []
        for motif in cached['motifs']:
            for parent in motif['parents']:
                data.append({
                    'ml_release_id': cached['release'],
                    'motif_id': motif['motif_id'],
                    'parent_ml_release_id': cached['parent'],
                    'parent_motif_id': parent['name']['full'],
                })
        return data

    def no_parents(self, data):
        """Check if the data has no parents. This is done through a different
        method than just attempting to create parents and seeing if the list is
        empty. This gives a method to determine if there are parents seperately
        from computing them for storage.

        Parameters
        ----------
        data : dict
            The cached grouping to examine

        Returns
        -------
        has_no_parents : bool
            True if there are no parents.
        """

        counts = data['parent_counts']['motifs']['unchanged'] + \
            data['parent_counts']['motifs']['updated']
        return not counts

    def data(self, pair, **kwargs):
        """Compute the parentage data. This will raise a skip exception if
        there are not parents, or if this is the first release (parent release
        is the same as the current release). This requires that there is data
        stored in the NR_CACHE_NAME file. If there is not, then an exception is
        raised.

        Parameters
        ----------
        release : str
            The nr release id to process.

        Raises
        ------
        Skip
            If this is the first release, or there are no parents.

        Returns
        -------
        data : list
            A list of dicts that can be written to the ml_parents table.
        """

        loop_type, release = pair
        cached = self.cached(loop_type)
        if not cached:
            raise core.InvalidState("No cached data")

        if cached['release'] == cached['parent']:
            raise core.Skip("No parents for first release")
        if self.no_parents(cached):
            raise core.Skip("Parent counts show no parents")

        return self.parents(cached)
