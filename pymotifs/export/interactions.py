from pymotifs import core
from pymotifs.models import UnitPairsInteractions
from pymotifs.interactions.pairwise import Loader as InterLoader
from pymotifs.utils import row2dict


class Exporter(core.Exporter):
    headers = ['unit_id1', 'unit_id2', 'FR3D basepair (f_lwbp)',
               'FR3D stacking (f_stacks)', 'FR3D base phosphate (f_bphs)']
    dependencies = set([InterLoader])
    compressed = True

    def filename(self, *args, **kwargs):
        return self.config['locations']['interactions_gz']

    def interactions(self, pdb):
        """Lookup all interactions for the given structure.
        """

        with self.session() as session:
            query = session.query(
                UnitPairsInteractions.unit_id_1.label(self.headers[0]),
                UnitPairsInteractions.unit_id_2.label(self.headers[1]),
                UnitPairsInteractions.f_lwbp.label(self.headers[2]),
                UnitPairsInteractions.f_stacks.label(self.headers[3]),
                UnitPairsInteractions.f_bphs.label(self.headers[4])
            ).filter_by(pdb_id=pdb)

        count = query.count()
        if not count:
            self.logger.warning("No interactions found for %s", pdb)
        else:
            self.logger.info("Found %s interactions for %s", count, pdb)

        return [row2dict(result) for result in query]

    def data(self, pdbs, **kwargs):
        """Load all interactions for the given structure. This returns a
        generator over all interactions.
        """

        for pdb in pdbs:
            self.logger.info('Writing out interactions for %s', pdb)
            for interaction in self.interactions(pdb):
                yield interaction
