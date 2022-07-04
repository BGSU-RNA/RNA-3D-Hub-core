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
        return session.query(mod.PairAnnotations).filter_by(pdb_id=pdb)
        #return session.query(mod.PairAnnotations).filter_by(pdb_id=pdb).filter_by(category='basepair')

    def to_process(self, pdbs, **kwargs):
        """
        Return PDB ids which do not have annotations in the table.

        :param list pdb: The list of pdb ids. Currently ignored.
        :param dict kwargs: The keyword arguments which are ignored.
        :returns: A list of correspondence ids to process.
        """

        with self.session() as session:
            query = session.query(mod.PairAnnotations.pdb_id).distinct()
            annotated = set([result.pdb_id for result in query])

        #pdb_list = sorted(list(set(pdbs) - annotated),reverse=True)
        pdb_list = sorted(list(set(pdbs) - annotated))

        # when adding a new category of annotations, you can process all PDBs
        #pdb_list = pdbs

        # if you need to run it twice, split the job in half
        #pdb_list = pdb_list[:4000]
        #pdb_list = pdb_list[4002:]

        if len(pdb_list) == 0:
            raise core.Skip("No PDB files need pair interactions annotated")

        return pdb_list


    def data(self, pdb, **kwargs):

        categories = {}
        categories['sO'] = []
        categories['basepair'] = []
        categories['coplanar'] = []
        categories['stacking'] = []

        structure = self.structure(pdb)
        interaction_to_triple_list, category_to_interactions, timerData, pair_to_data = annotate_nt_nt_in_structure(structure,categories)

        annotations = []

        all_categories = list(set(categories.keys()) & set(category_to_interactions.keys()))

        for category in all_categories:
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
                        'unit_id_1' : triple[0],
                        'unit_id_2' : triple[1],
                        'annotation': interaction_simple,
                        'annotation_detail' : interaction_detail,
                        'category'  : category,
                        'crossing'  : triple[2]
                        })

        if len(annotations) == 0:
            annotations.append({
                'pdb_id'    : pdb,
                'unit_id_1' : 'placeholder',
                'unit_id_2' : 'placeholder',
                'annotation': None,
                'annotation_detail' : None,
                'category'  : 'placeholder',
                'crossing'  : None
                })

        return annotations
