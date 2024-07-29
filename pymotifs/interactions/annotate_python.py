"""
Run Python code to annotate pairwise interactions
and store them in unit_pairs_annotations, noting that program = 'python'
"""

from pymotifs import core
from pymotifs import models as mod

from pymotifs.pdbs.info import Loader as PdbLoader
from pymotifs.units.info import Loader as UnitLoader

# import these from the interactions directory; in the future, they'll be updated with fr3d-python
from fr3d.classifiers.NA_pairwise_interactions import annotate_nt_nt_in_structure

class Loader(core.SimpleLoader):
    allow_no_data = True
    dependencies = set([PdbLoader, UnitLoader])

    @property
    def table(self):
        return mod.UnitPairsInteractions

    def query(self, session, pdb):
        return session.query(mod.UnitPairsInteractions).\
                filter_by(pdb_id=pdb).\
                filter_by(program='python')

    def to_process(self, pdbs, **kwargs):
        """
        Return PDB ids which do not have python annotations in the table.

        Until we want to annotate every RNA structure with Python annotations,
        exclude those structures with Matlab annotations.

        :param list pdbs: The list of pdb ids of interest on this run.
        :param dict kwargs: The keyword arguments which are ignored.
        :returns: A list of pdb ids to process.
        """

        if len(pdbs) < 100:
            # get (pdb_id,program) pairs in unit_pairs_interactions for these pdbs only; faster
            with self.session() as session:
                query = session.query(mod.UnitPairsInteractions.pdb_id,mod.UnitPairsInteractions.program).\
                filter(mod.UnitPairsInteractions.pdb_id.in_(pdbs)).\
                distinct()
                pdb_id_program_all = set([(result.pdb_id,result.program) for result in query])
        else:
            # get all (pdb_id,program) pairs in unit_pairs_interactions table
            # sometimes fast, sometimes not
            with self.session() as session:
                query = session.query(mod.UnitPairsInteractions.pdb_id,mod.UnitPairsInteractions.program).\
                distinct()
                pdb_id_program_all = set([(result.pdb_id,result.program) for result in query])

        self.logger.info('Found %d pdb_id,program pairs' % len(pdb_id_program_all))

        pdbs_matlab = set([pdb_id for pdb_id,program in pdb_id_program_all if program == 'matlab'])
        pdbs_python = set([pdb_id for pdb_id,program in pdb_id_program_all if program == 'python'])

        need_to_annotate = set(pdbs) - pdbs_python - pdbs_matlab

        if len(need_to_annotate) == 0:
            raise core.Skip("No PDB files need pairwise interactions annotated with python")
        else:
            return sorted(need_to_annotate)


    def data(self, pdb, **kwargs):
        """
        read the .cif file for the pdb id, annotate interactions, store them
        """

        # interaction categories to process; names from NA_pairwise_interactions.py
        # basepair,basepair_detail,stacking,so,backbone,coplanar,covalent,sugar_ribose,near,bss,loops
        categories = {}
        categories['basepair'] = []
        categories['basepair_detail'] = []
        categories['stacking'] = []
        categories['sO'] = []
        categories['backbone'] = []
        categories['coplanar'] = []
        categories['covalent'] = []
        categories['sugar_ribose'] = []
        categories['bss'] = []
        categories['loops'] = []

        category_to_heading = {}
        category_to_heading['basepair'] = 'f_lwbp'
        category_to_heading['basepair_detail'] = 'f_lwbp_detail'
        category_to_heading['stacking'] = 'f_stacks'
        category_to_heading['backbone'] = ''  # need special processing for BPh, BR
        category_to_heading['sO'] = 'f_so'
        category_to_heading['coplanar'] = 'f_coplanar'
        category_to_heading['sugar_ribose'] = 'f_sugar_ribose'
        category_to_heading['bss'] = 'f_bss'
        category_to_heading['covalent'] = 'f_covalent'

        # load the .cif file
        structure = self.structure(pdb)

        # annotate interactions
        interaction_to_triple_list, category_to_interactions, timerData, pair_to_data = annotate_nt_nt_in_structure(structure,categories)

        # actual categories found
        all_categories = list(set(categories.keys()) & set(category_to_interactions.keys()))

        pair_to_dictionary = {}

        for category in all_categories:
            heading = category_to_heading[category]

            for interaction in category_to_interactions[category]:
                # if category == 'basepair':
                #     if len(interaction) == 3:
                #         interaction_simple = interaction[0] + interaction[1:3].upper()
                #     elif len(interaction) == 4:
                #         interaction_simple = interaction[0:2] + interaction[2:4].upper()
                #     else:
                #         interaction_simple = interaction
                #     interaction_detail = interaction
                # else:
                #     interaction_simple = interaction
                #     interaction_detail = None



                if category == 'backbone':
                    if 'BPh' in interaction:
                        heading = 'f_bphs'
                    elif 'BR' in interaction:
                        heading = 'f_brbs'

                for triple in interaction_to_triple_list[interaction]:
                    u1 = triple[0]
                    u2 = triple[1]
                    pair = u1 + "+" + u2

                    # self.logger.info('Adding interaction %s between %s and %s in %s' % (interaction,u1,u2,pdb))

                    if not pair in pair_to_dictionary:
                        pair_to_dictionary[pair] = {}

                    pair_to_dictionary[pair]['unit_id_1'] = u1
                    pair_to_dictionary[pair]['unit_id_2'] = u2
                    pair_to_dictionary[pair]['pdb_id'] = pdb
                    pair_to_dictionary[pair]['program'] = 'python'
                    pair_to_dictionary[pair][heading] = interaction

                    if 'f_crossing' in pair_to_dictionary[pair]:
                        if not pair_to_dictionary[pair]['f_crossing'] == triple[2]:
                            print('Different crossing numbers for %s and %s interaction %s' % (u1,u2,interaction))
                    else:
                        pair_to_dictionary[pair]['f_crossing'] = triple[2]


        # list of annotations to return
        annotations = []
        for pair in sorted(pair_to_dictionary.keys()):
            annotations.append(pair_to_dictionary[pair])

        """
        # if we want to leave a note that we checked for interactions, we could try this
        if len(pair_to_dictionary) == 0:
            annotations.append({
                'pdb_id'    : pdb,
                'unit_id_1' : 'placeholder',
                'unit_id_2' : 'placeholder',
                'annotation': None,
                'annotation_detail' : None,
                'category'  : 'placeholder',
                'crossing'  : None
                })
        """

        if len(annotations) == 0:
            raise core.Skip("No new interactions to add for %s" % pdb)

        return annotations
