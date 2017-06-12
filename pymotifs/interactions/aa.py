"""Import NA/AA interactions
"""

import operator as op

from pymotifs import core
from pymotifs import models as mod

from pymotifs.units.info import Loader as UnitLoader
from pymotifs.pdbs.info import Loader as PdbLoader

from fr3d.classifiers.base_aafg import Classifier


class Loader(core.SimpleLoader):
    allow_no_data = True
    dependencies = set([UnitLoader, PdbLoader])

    @property
    def table(self):
        return mod.UnitAaInteractions

    def query(self, session, pdb):
        return session.query(self.table).filter_by(pdb_id=pdb)

    def annotations(self, structure):
        classifier = Classifier()
        pairs = []
        for na_id, aa_id, (annotation, value) in classifier.classify(structure):
            pairs.append({
                'na_unit_id': na_id,
                'aa_unit_id': aa_id,
                'pdb_id': structure.pdb,
                'annotation': annotation,
                'value': value,
            })
        return sorted(pairs, key=op.itemgetter('na_unit_id', 'aa_unit_id'))

    def data(self, pdb, **kwargs):
        structure = self.structure(pdb)
        structure.infer_hydrogens()
        return self.annotations(structure)
