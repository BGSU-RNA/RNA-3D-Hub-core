"""Stage to populate the units.info table.
The stage only need to run once
"""

import itertools as it

import pymotifs.core as core

from pymotifs import models as mod
from pymotifs.utils import units
from pymotifs import utils
# from pymotifs.download import Downloader
# from pymotifs.pdbs.info import Loader as PdbLoader
# from pymotifs.units import Loader as UnitInfoLoder
from sqlalchemy import and_
from collections import defaultdict
from Bio.Alphabet import ThreeLetterProtein

AA = [seq.upper() for seq in ThreeLetterProtein().letters]

class Loader(core.Loader):
    merge_data = True
    mark = False
    allow_no_data = True


    dependencies = set([])
    """The dependencies for this stage."""

    def to_process(self, pdbs, **kwargs):
        with self.session() as session:
            query = session.query(mod.UnitRotations.unit_id,mod.UnitRotations.pdb_id).\
                    filter(mod.UnitRotations.pdb_id.in_(['1ASY', '1EVV']))
        d = defaultdict(dict)
        for result in query:
            if result.unit_id.split('|')[3] not in ['A', 'C', 'G', 'U','DA', 'DC', 'DG', 'DT']:
                if d.get(result.pdb_id):
                    d[result.pdb_id].append(result.unit_id)
                else:
                    d[result.pdb_id] = []
                    d[result.pdb_id].append(result.unit_id)
        if not len(d):
                raise core.Skip("Skipping summary, no new correspondences")
        return list(d.values())

    def remove(self, pdb_id, **kwargs):
        self.logger.info("Not removing anything")

    def has_data(self, unit_id, **kwargs):
        return 0


    def data(self, unit_ids, **kwargs):
        """Get the rotation matrices for all RNA residues in the given pdb.

        :pdb: The pdb to process
        :yields: Yields a series of rotation matrices.
        """
        structure = self.structure(unit_ids[0][:4])
        for residue in structure.residues():
            if residue.unit_id() in unit_ids:
                if hasattr(residue, 'rotation_matrix'):
                    matrix = residue.rotation_matrix
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
                    # print(data)
                    yield mod.UnitRotations(**data)

