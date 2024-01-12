"""
Run Python code to annotate pairwise interactions
and store them in unit_pairs_annotations alongside
the Matlab annotations.
"""

import operator as op

from pymotifs import core
from pymotifs import models as mod

from pymotifs.units.info import Loader as UnitLoader
from pymotifs.pdbs.info import Loader as PdbLoader

# import these from the interactions directory; in the future, they'll be updated with fr3d-python
from fr3d.classifiers.NA_pairwise_interactions import annotate_nt_nt_in_structure
from fr3d.classifiers.class_limits_2023 import nt_nt_cutoffs   # use latest cutoffs

class Loader(core.SimpleLoader):
    allow_no_data = True
    dependencies = set([UnitLoader, PdbLoader])

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

        :param list pdbs: The list of pdb ids of interest on this run.
        :param dict kwargs: The keyword arguments which are ignored.
        :returns: A list of pdb ids to process.
        """

        self.logger.info('Starting to_process query')

        # structures that are already annotated by python
        with self.session() as session:
            query = session.query(mod.UnitPairsInteractions.pdb_id).\
            filter(mod.UnitPairsInteractions.program == 'python').\
            distinct()

            annotated = set([result.pdb_id for result in query])

        self.logger.info('Found %d pdb ids already annotated by python' % len(annotated))

        # temporary:  DNA-containing structures
        with self.session() as session:
            query = session.query(mod.ChainInfo.pdb_id).\
            filter(mod.ChainInfo.macromolecule_type == 'DNA').\
            distinct()

            DNA = set([result.pdb_id for result in query])

        DNA_list = sorted(DNA - annotated)

        return DNA_list



        pdb_list = sorted(list(set(pdbs) - annotated))
        #pdb_list = sorted(list(set(pdbs) - annotated), reverse = True)

        # when adding a new category of annotations, you can process all PDBs
        #pdb_list = pdbs

        # if you need to run it twice, split the job roughly in half and run in different directions
        #pdb_list = pdb_list[:3800]
        #pdb_list = pdb_list[3800:]

        if len(pdb_list) == 0:
            raise core.Skip("No PDB files need pair interactions annotated")

        return pdb_list


    def data(self, pdb, **kwargs):
        """
        read the .cif file for the pdb id, annotate interactions, store
        """

        # future:  query the unit_pairs_interactions table to find interactions with this pdb id that are already stored
        # future:  and only record new interactions

        # interaction categories to process
        categories = {}
        categories['sO'] = []
        categories['basepair'] = []
        categories['coplanar'] = []
        categories['stacking'] = []

        category_to_heading = {}
        category_to_heading['basepair'] = 'f_lwbp'
        category_to_heading['stacking'] = 'f_stacks'
        category_to_heading['backbone'] = 'backbone'
        category_to_heading['sO'] = 'f_so'
        category_to_heading['coplanar'] = 'f_coplanar'

        # load the .cif file
        structure = self.structure(pdb)

        # annotate
        interaction_to_triple_list, category_to_interactions, timerData, pair_to_data = annotate_nt_nt_in_structure(structure,categories)

        all_categories = list(set(categories.keys()) & set(category_to_interactions.keys()))

        pair_to_dictionary = {}

        for category in all_categories:
            heading = category_to_heading[category]

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

                if 'BPh' in interaction:
                    heading = 'f_bphs'
                elif 'BR' in interaction:
                    heading = 'f_brbs'

                for triple in interaction_to_triple_list[interaction]:
                    u1 = triple[0]
                    u2 = triple[1]
                    pair = u1 + "+" + u2

                    f1 = u1.split("|")
                    f2 = u2.split("|")

                    b1 = f1[3]
                    b2 = f2[3]

                    # temporary in 2023-02-27: skip anything Matlab may have already annotated
                    if b1 in ['A','C','G','U'] and b2 in ['A','C','G','U']:
                        continue

                    if not pair in pair_to_dictionary:
                        pair_to_dictionary[pair] = {}

                    pair_to_dictionary[pair]['unit_id_1'] = u1
                    pair_to_dictionary[pair]['unit_id_2'] = u2
                    pair_to_dictionary[pair]['pdb_id'] = pdb
                    pair_to_dictionary[pair][heading] = interaction_simple

                    if 'f_crossing' in pair_to_dictionary[pair]:
                        if not pair_to_dictionary[pair]['f_crossing'] == triple[2]:
                            print('Different crossing numbers for %s and %s interaction %s' % (u1,u2,interaction_simple))
                    else:
                        pair_to_dictionary[pair]['f_crossing'] = triple[2]

                    pair_to_dictionary[pair]['program'] = 'python'

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
