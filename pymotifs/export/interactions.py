"""
Export the interactions for all structures to a single .gz file.
If given structures it will only export interactions from those
pdbs but if given all it will export all interactions.
"""

import itertools as it

from pymotifs import core
from pymotifs import models as mod
from pymotifs.interactions.pairwise import Loader as InterLoader
from pymotifs.utils import row2dict


class Exporter(core.Exporter):

    headers = ['unit_id1', 'unit_id2', 'FR3D basepair (f_lwbp)',
               'FR3D stacking (f_stacks)', 'FR3D base phosphate (f_bphs)']
    dependencies = set([InterLoader])
    compressed = True
    mark = False

    def filename(self, *args, **kwargs):
        """This will always return the configured path at
        locations.interactions_gz and is where the interaction export will be
        written.

        Returns
        -------
        filename : str
            The path to write to.
        """
        return self.config['locations']['interactions_gz']

    def to_process(self, pdbs, **kwargs):
        """

        """

        if len(pdbs) < 500:
            raise core.Skip("Too few pdb files being processed to write all interactions")

        return pdbs

    def interactions(self, pdb):
        """Look up all interactions for the given structure. This gets all
        interaction entries. If there are none this returns an empty list. The
        entries in the list are dictonaries with the same names as in
        `Exporter.headers`.

        Parameters
        ----------
        pdb : str
            The PDB id to look up interactions for

        Returns
        -------
        interactions : list
            A list of all interactions.
        """

        with self.session() as session:
            query = session.query(
                mod.UnitPairsInteractions.unit_id_1.label(self.headers[0]),
                mod.UnitPairsInteractions.unit_id_2.label(self.headers[1]),
                mod.UnitPairsInteractions.f_lwbp.label(self.headers[2]),
                mod.UnitPairsInteractions.f_stacks.label(self.headers[3]),
                mod.UnitPairsInteractions.f_bphs.label(self.headers[4])
            ).filter_by(pdb_id=pdb,program='matlab')

            count = query.count()
            self.logger.info("Found %5d interactions for %s" % (count, pdb))

            return [row2dict(result) for result in query]

    def data(self, pdbs, **kwargs):
        """Load all interactions for the given structure. This returns a
        generator over all interactions.

        Parameters
        ----------
        pdbs : list
            The pdbs to look up.

        Returns
        -------
        interactions : iterable
            An iterable over all interactions in all the given structures.
        """

        if len(pdbs) < 500:
            raise core.Skip('Too few pdb files being processed to export interactions')

        interactions = it.imap(self.interactions, pdbs)
        interactions = it.chain.from_iterable(interactions)
        return interactions
