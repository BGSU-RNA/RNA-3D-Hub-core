"""
Find units that do not have a glycosidic center and/or rotation matrix calculated, and load them.
Restrict to modified nucleotides if you want.

bin/pipeline.py --log-file logs/units.rotation_by_unit_id.log --log-mode w run --skip-dependencies --all units.rotation_by_unit_id
"""

import glob
import os

import pymotifs.core as core
from pymotifs import models as mod
from pymotifs.units.info import Loader as InfoLoader

class Loader(core.SimpleLoader):
    """
    A class to load missing centers and rotation matrices into the database.
    """

    dependencies = set([InfoLoader])
    allow_no_data = True
    mark = False

    modified_only = False
    modified_only = True

    def to_process(self, pdbs, **kwargs):
        """
        If pdbs is a list of just one pdb id, process that one.
        Otherwise, find all pdbs that do not have rotation matrix
        """

        if len(pdbs) == 1:
            # ok, this does not quite work.  Fix it if you need it.

            # find all unit ids that are marked as rna or dna
            # when new modified nucleotides are added, we'll need to add rna/dna
            with self.session() as session:
                query = session.query(mod.UnitInfo.unit_id).\
                                filter(mod.UnitInfo.unit_type_id.in_(['rna','dna'])).\
                                filter(mod.UnitInfo.pdb_id == pdbs[0])

                na_units = set([r.unit_id for r in query])
                self.logger.info('Found %d na units' % len(na_units))

            # find all unit ids with rotation matrices
            with self.session() as session:
                query = session.query(mod.UnitRotations.unit_id).\
                                filter(mod.UnitInfo.pdb_id == pdbs[0])

                rotation_units = set([r.unit_id for r in query])
                self.logger.info('Found %d units with rotation matrices' % len(rotation_units))

        else:
            # find all unit ids that are marked as rna or dna
            # when new modified nucleotides are added, we'll need to add rna/dna
            # you could list the new modified nucleotides here, to target them
            if self.modified_only:
                # exclude standard bases and common nucleotides with no base atoms
                # with self.session() as session:
                #     query = session.query(mod.UnitInfo.unit_id).\
                #                     filter(mod.UnitInfo.unit_type_id.in_(['rna','dna'])).\
                #                     filter(~mod.UnitInfo.unit.in_(['A','C','G','U','DA','DC','DG','DT','N','3DR']))
                #     na_units = set([r.unit_id for r in query])
                #     self.logger.info('Found %d modified na units' % len(na_units))

                # target specific new modified nucleotides
                with self.session() as session:
                    query = session.query(mod.UnitInfo.unit_id).\
                                    filter(mod.UnitInfo.unit_type_id.in_(['rna','dna'])).\
                                    filter(mod.UnitInfo.unit.in_(['DU','DI','A1EFN','B4P','A1MA9','A1L89','A1L82','A1L3P','A1B8A','ZDU','USM']))
                    na_units = set([r.unit_id for r in query])
                    self.logger.info('Found %d modified na units' % len(na_units))
            else:
                # query takes about 90 seconds
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

            # find all unit ids with rotation matrices
            with self.session() as session:
                query = session.query(mod.UnitRotations.unit_id)

                rotation_units = set([r.unit_id for r in query])
                self.logger.info('Found %d units with rotation matrices' % len(rotation_units))

        units_without_glycosidic = na_units - glycosidic_units
        units_without_rotation = na_units - rotation_units

        self.logger.info('Found %d units lacking glycosidic center' % len(units_without_glycosidic))
        self.logger.info('Found %d units lacking rotation matrix  ' % len(units_without_rotation))

        if len(units_without_glycosidic) == 0 and len(units_without_rotation) == 0:
            raise core.Skip("All na units have glycosidic center and rotation matrix")

        # find all pdb ids with at least one unit id needing work
        pdb_to_units = {}
        for u in units_without_glycosidic:
            fields = u.split("|")
            pdb = fields[0]
            if not pdb in pdb_to_units:
                pdb_to_units[pdb] = {'glycosidic': [], 'rotation': []}
            pdb_to_units[pdb]['glycosidic'].append(u)
        for u in units_without_rotation:
            fields = u.split("|")
            pdb = fields[0]
            if not pdb in pdb_to_units:
                pdb_to_units[pdb] = {'glycosidic': [], 'rotation': []}
            pdb_to_units[pdb]['rotation'].append(u)

        if len(pdb_to_units) == 0:
            raise core.Skip("No pdbs need glycosidic center or rotation matrices added")

        self.logger.info('Found %d pdb files with missing glycosidic center and/or rotation matrix' % len(pdb_to_units))

        # return (pdb,units) tuples to process, starting with most recent PDB files
        return sorted(pdb_to_units.items(),reverse=True)


    def query(self, session, pdb_units):
        # always return an empty query, because we already
        # focus on unit ids that need to be processed
        return session.query(mod.UnitRotations).\
            filter(mod.UnitRotations.pdb_id == 'XX')


    # def type(self, unit):
    #     """Compute the component type, ie A, C, G, U is RNA, DA, DC, etc is DNA
    #     and so forth.

    #     Parameters
    #     ----------
    #     unit : Component
    #         The unit to get the component for

    #     Returns
    #     -------
    #     component_type : str
    #         The component type.
    #     """
    #     return units.component_type(unit)


    def data(self, pdb_units, **kwargs):

        pdb = pdb_units[0]
        glycosidic_needed = pdb_units[1]['glycosidic']
        rotation_needed = pdb_units[1]['rotation']

        # load the structure and get the residues
        structure = self.structure(pdb)
        for residue in structure.residues():
            if residue.unit_id() in glycosidic_needed:
                for name in residue.centers.definitions():
                    center = residue.centers[name]
                    if len(center) == 3 and name == 'glycosidic':
                        self.logger.info('Adding glycosidic center for %s' % residue.unit_id())
                        # delete pickle files with centers and rotations,
                        # so they will be exported again the next time
                        for f in glob.glob("/usr/local/pipeline/hub-core/data/units/%s-*NA.pickle" % pdb):
                            os.remove(f)
                        yield mod.UnitCenters(unit_id=residue.unit_id(),
                                            name=name,
                                            pdb_id=pdb,
                                            x=float(center[0]),
                                            y=float(center[1]),
                                            z=float(center[2]))

            if residue.unit_id() in rotation_needed:
                print('Need rotation matrix for %s' % residue.unit_id())
                self.logger.info('Need rotation matrix for %s' % residue.unit_id())
                if hasattr(residue, 'rotation_matrix'):
                    matrix = residue.rotation_matrix
                    self.logger.info(str(matrix))
                    if matrix is not None and matrix[0, 0] is not None:
                        self.logger.info('Adding rotation matrix for %s' % residue.unit_id())
                        # delete pickle files with centers and rotations,
                        # so they will be exported again the next time
                        for f in glob.glob("/usr/local/pipeline/hub-core/data/units/%s-*NA.pickle" % pdb):
                            os.remove(f)
                        data = {'unit_id':residue.unit_id(),
                                'pdb_id':residue.unit_id()[:4],
                                'cell_0_0':float(matrix[0, 0]),
                                'cell_0_1':float(matrix[0, 1]),
                                'cell_0_2':float(matrix[0, 2]),
                                'cell_1_0':float(matrix[1, 0]),
                                'cell_1_1':float(matrix[1, 1]),
                                'cell_1_2':float(matrix[1, 2]),
                                'cell_2_0':float(matrix[2, 0]),
                                'cell_2_1':float(matrix[2, 1]),
                                'cell_2_2':float(matrix[2, 2])}
                        yield mod.UnitRotations(**data)
                    else:
                        self.logger.info('No rotation matrix for %s' % residue.unit_id())
                else:
                    self.logger.info('No rotation matrix for %s' % residue.unit_id())
