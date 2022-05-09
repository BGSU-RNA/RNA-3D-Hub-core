"""
Run Python code to annotate pairwise interactions
"""

import operator as op

from pymotifs import core
from pymotifs import models as mod

from pymotifs.units.info import Loader as UnitLoader
from pymotifs.pdbs.info import Loader as PdbLoader

from NA_pairwise_interactions import annotate_nt_nt_in_structure

class Loader(core.SimpleLoader):
    allow_no_data = True
    dependencies = set([UnitLoader, PdbLoader])

    @property
    def table(self):
        return mod.PairAnnotations

    def query(self, session, pdb):
        return session.query(self.table).filter_by(pdb_id=pdb)

    def must_recompute(self, *args, **kwargs):
        return True

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

        categories = {}
        categories['sO'] = 'sO'

        structure = self.structure(pdb)
        interaction_to_triple_list, category_to_interactions = annotate_nt_nt_in_structure(structure,categories)

        annotation_types = categories.keys()

        annotations = []

        for category in categories.keys():
            for interaction in category_to_interactions[category]:
                if category == 'basepair':
                    if len(interaction) == 3:
                        interaction_simple = interaction[0] + interaction[1:3].upper()
                    elif len(interaction) == 4:
                        interaction_simple = interaction[0:2] + interaction[2:4].upper()
                    else:
                        interaction_simple = interaction
                    interaction_detail = interaction
                else:
                    interaction_simple = interaction
                    interaction_detail = None

                for triple in interaction_to_triple_list[interaction]:
                    annotations.append({
                        'pdb_id'    : pdb,
                        'unit_id_1' : triple[1],
                        'unit_id_2' : triple[2],
                        'annotation': interaction_simple,
                        'annotation_detail' : interaction_detail,
                        'category'  : category,
                        'crossing'  : triple[3]
                        })

        return annotations
