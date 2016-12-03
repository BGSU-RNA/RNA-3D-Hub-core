"""Write the loop to motif assignment data. This will populate the ml_loops
table in the database. This will load the cached data to store all motifs into
the DB.
"""

from pymotifs import core
from pymotifs import models as mod

from pymotifs.motifs.utils import BaseLoader
from pymotifs.motifs.info import Loader as InfoLoader
from pymotifs.loops.extractor import Loader as LoopLoader
from pymotifs.motifs.release import Loader as ReleaseLoader


class Loader(BaseLoader):
    dependencies = set([ReleaseLoader, LoopLoader, InfoLoader])

    @property
    def table(self):
        return mod.MlLoops

    def assignments(self, cached):
        """Compute the assignments of loop to motif.

        Parameters
        ----------
        cached : dict
            The cached motif data.

        Returns
        -------
        assignments : list
            A list of loop to motif assignments in the form of MlLoops objects.
        """
        data = []
        for entry in cached['motifs']:
            for loop in entry['ordering']:
                data.append(self.table(
                    loop_id=loop['loop_id'],
                    motif_id=entry['motif_id'],
                    ml_release_id=cached['release'],
                ))
        return data

    def data(self, pair, **kwargs):
        """Compute the loop to motif data to store.

        Parameters
        ----------
        pair : (loop_type, release)
            The tuple of the loop type and release to load.

        Returns
        -------
        assignments : list
            A list of loop to motif assignments in the form of MlLoops objects.
        """
        loop_type, release = pair
        cached = self.cached(loop_type)
        if not cached:
            raise core.InvalidState("No cached data")
        return self.assignments(cached)
