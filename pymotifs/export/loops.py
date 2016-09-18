"""A loader to write the loop export file. This will create a CSV file that is
provided to NDB for information about all loops in all structures.
"""

import itertools as it

from pymotifs import core
from pymotifs import models as mod

from pymotifs.utils import row2dict

from pymotifs.loops.extractor import Loader as LoopLoader
from pymotifs.loops.positions import Loader as PositionLoader


class Exporter(core.Exporter):
    """The actual stage that gets run."""

    headers = ['id', 'motif_id', 'pdb', 'nts']
    dependencies = set([LoopLoader, PositionLoader])

    compressed = True
    """We provide a compressed file."""

    mark = False

    def filename(self, *args, **kwargs):
        """The filename to write to. It is always the configured
        locations.loop_gz file.
        """
        return self.config['locations']['loops_gz']

    def current_ml_release(self):
        """Fetch the current ml release. If there is no ml_release_id then we
        will return 0.0.

        Returns
        -------
        ml_release : str
            The current ml_release_id.
        """

        with self.session() as session:
            current = session.query(mod.MlReleases.ml_releases_id).\
                order_by(mod.MlReleases.date).\
                limit(1).\
                first()

            if not current:
                return '0.0'

            return current.ml_release_id

    def loops(self, pdb):
        """Get all loops in the current structure. If the loop is part of the
        current motif atlas release we will fetch the motif assignment as well.

        Parameters
        ----------
        pdb : str
            The pdb id to look up structures for.

        Returns
        -------
        loops : list
            A list of loop dictonaries that contain an 'id', 'pdb', 'nts' and
            'motif_id' column.
        """

        current_ml_release = self.current_ml_release()
        with self.session() as session:
            query = session.query(mod.LoopInfo.loop_id.label('id'),
                                  mod.LoopInfo.pdb_id.label('pdb'),
                                  mod.LoopInfo.unit_ids.label('nts'),
                                  mod.MlLoops.motif_id.label('motif_id')
                                  ).\
                outerjoin(mod.MlLoops,
                          (mod.MlLoops.loop_id == mod.LoopInfo.loop_id) &
                          (mod.MlLoops.ml_release_id == current_ml_release)).\
                filter(mod.LoopInfo.pdb_id == pdb).\
                order_by(mod.LoopInfo.loop_id)

            count = query.count()
            if not count:
                self.logger.info("No loops found for %s", pdb)
            else:
                self.logger.info("Found %s loops for %s", count, pdb)

            return [row2dict(result) for result in query]

    def data(self, pdbs, **kwargs):
        """Get the loop data for all structures.

        Parameters
        ----------
        pdbs : list
            The pdb ids to look up.

        Returns
        -------
        loops : iterator
            An iterator over loop dictonary as from `Exporter.loop`.
        """

        loops = it.imap(self.loops, pdbs)
        loops = it.chain.from_iterable(loops)
        return loops
