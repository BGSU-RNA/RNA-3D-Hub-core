"""
Find units that do not have glycosidic centers calculated, and load them.
bin/pipeline.py --log-file logs/units.centers_glycosidic_by_unit_id.log --log-mode w run --skip-dependencies --all units.centers_glycosidic_by_unit_id

"""

import pymotifs.core as core
from pymotifs import models as mod
from pymotifs.units.info import Loader as InfoLoader
from pymotifs.skip_files import SKIP
from pymotifs.utils import units

class Loader(core.SimpleLoader):
    """
    A class to load all glycosidic centers into the database.
    """

    dependencies = set([InfoLoader])
    allow_no_data = True
    mark = False

    def to_process(self, pdbs, **kwargs):
        """
        If pdbs is a list of just one pdb id, process that one.
        Otherwise, find all rna, dna units that do not have glycosidic
        """

        if len(pdbs) == 1:
            return pdbs
        else:
            # find all unit ids that are marked as rna or dna
            # when new modified nucleotides are added, we'll need to add rna/dna
            # takes about 90 seconds
            with self.session() as session:
                query = session.query(mod.UnitInfo.unit_id).\
                               filter(mod.UnitInfo.unit_type_id.in_(['rna','dna']))

                na_units = set([r.unit_id for r in query])
                self.logger.info('Found %d na units' % len(na_units))

            # find all unit ids with glycosidic centers
            # takes almost 3 minutes
            with self.session() as session:
                query = session.query(mod.UnitCenters.unit_id).\
                               filter(mod.UnitCenters.name == 'glycosidic')

                glycosidic_units = set([r.unit_id for r in query])
                self.logger.info('Found %d units with glycosidic centers' % len(glycosidic_units))

            units_without_glycosidic = na_units - glycosidic_units

            if len(units_without_glycosidic) == 0:
                raise core.Skip("All na units have glycosidic centers")

            # find all pdb ids with at least one unit id in units_without_glycosidic
            pdbs_to_load = set([u.split("|")[0] for u in units_without_glycosidic]) - set(SKIP)

            if len(pdbs_to_load) == 0:
                raise core.Skip("No pdbs need glycosidic centers added")

            return sorted(pdbs_to_load)


    def query(self, session, pdb):
        return session.query(mod.UnitCenters).\
            filter(mod.UnitCenters.pdb_id == pdb).\
            filter(mod.UnitCenters.name == 'glycosidic')


    def type(self, unit):
        """Compute the component type, ie A, C, G, U is RNA, DA, DC, etc is DNA
        and so forth.

        Parameters
        ----------
        unit : Component
            The unit to get the component for

        Returns
        -------
        component_type : str
            The component type.
        """
        return units.component_type(unit)


    def data(self, pdb, **kwargs):

        # load the structure
        structure = self.structure(pdb)
        for residue in structure.residues():
            for name in residue.centers.definitions():
                center = residue.centers[name]
                if len(center) == 3 and name == 'glycosidic':
                    self.logger.info('Adding glycosidic center for %s' % residue.unit_id())
                    yield mod.UnitCenters(unit_id=residue.unit_id(),
                                          name=name,
                                          pdb_id=pdb,
                                          x=float(center[0]),
                                          y=float(center[1]),
                                          z=float(center[2]))
