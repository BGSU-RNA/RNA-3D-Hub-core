"""
Run Python code to annotate pairwise interactions
and store them in unit_pairs_annotations_2024, noting that program = 'fr3d'
2024 in the table name refers to the year when the parameters were set
Making a separate table will hopefully speed up queries for pairwise interactions
Next time there is a big change to the annotation program, make a new table
"""

import os

from pymotifs import core
from pymotifs import models as mod

from pymotifs.pdbs.info import Loader as PdbLoader
from pymotifs.units.info import Loader as UnitLoader
from pymotifs.chains.info import Loader as ChainLoader
from pymotifs.constants import DATA_FILE_DIRECTORY

# import these from the interactions directory; in the future, they'll be updated with fr3d-python
from fr3d.classifiers.NA_pairwise_interactions import annotate_nt_nt_in_structure
from fr3d.classifiers.NA_unit_annotation import annotate_self_base_backbone

class Loader(core.SimpleLoader):
    allow_no_data = True
    dependencies = set([PdbLoader, UnitLoader, ChainLoader])

    @property
    def table(self):
        return mod.UnitPairsInteractions2024

    def query(self, session, pdb):
        return session.query(mod.UnitPairsInteractions2024).\
                filter_by(pdb_id=pdb).\
                filter_by(program='fr3d')

    def to_process(self, pdbs, **kwargs):
        """
        Return PDB ids which do not have python annotations in the table.
        Because of the nature of the annotations, even a file with only
        a few nucleotides should have at least one pairwise interaction.

        :param list pdbs: The list of pdb ids of interest on this run.
        :param dict kwargs: The keyword arguments which are ignored.
        :returns: A list of pdb ids to process.
        """

        if len(pdbs) < 100:
            # get (pdb_id,program) pairs in unit_pairs_interactions_2024 for these pdbs only; faster
            with self.session() as session:
                query = session.query(mod.UnitPairsInteractions2024.pdb_id,mod.UnitPairsInteractions2024.program).\
                filter(mod.UnitPairsInteractions2024.pdb_id.in_(pdbs)).\
                distinct()
                pdb_id_program_all = set([(result.pdb_id,result.program) for result in query])
        else:
            # get all (pdb_id,program) pairs in unit_pairs_interactions_2024 table
            with self.session() as session:
                query = session.query(mod.UnitPairsInteractions2024.pdb_id,mod.UnitPairsInteractions2024.program).\
                distinct()
                pdb_id_program_all = set([(result.pdb_id,result.program) for result in query])

        if False:
            # temporary to get all pdb files in unit_info (unit ids will exist and not cause foreign key problems)
            with self.session() as session:
                query = session.query(mod.UnitInfo.pdb_id).\
                distinct()
                pdbs = [result.pdb_id for result in query]
            print('Found %d pdb files in unit_info' % len(pdbs))

            # get a list of .cif files in /usr/local/pipeline/hub-core/
            cif_files = os.listdir('/usr/local/pipeline/hub-core/cif_files/')
            cif_files = [f for f in cif_files if f.endswith('.cif.gz')]
            cif_files = [f.replace('.cif.gz','') for f in cif_files]

        self.logger.info('Found %d pdb_id,program pairs' % len(pdb_id_program_all))

        pdbs_fr3d = set([pdb_id for pdb_id,program in pdb_id_program_all if program == 'fr3d'])

        need_to_annotate = set(pdbs) - pdbs_fr3d

        print('Found %d files to annotate' % len(need_to_annotate))

        if len(need_to_annotate) == 0:
            raise core.Skip("No PDB files need pairwise interactions annotated with fr3d-python")
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
        categories['so'] = []
        categories['backbone'] = []
        categories['coplanar'] = []
        categories['covalent'] = []
        categories['sugar_ribose'] = []
        categories['bss'] = []
        categories['loops'] = []             # these are stored in different tables, but get them now

        category_to_heading = {}
        category_to_heading['basepair'] = 'f_lwbp'
        category_to_heading['basepair_detail'] = 'f_lwbp_detail'
        category_to_heading['stacking'] = 'f_stacks'
        category_to_heading['backbone'] = ''  # need special processing for BPh, BR
        category_to_heading['sO'] = 'f_so'
        category_to_heading['so'] = 'f_so'    # we may have missed all the sO interactions because of sO vs so!
        category_to_heading['coplanar'] = 'f_coplanar'
        category_to_heading['sugar_ribose'] = 'f_sugar_ribose'
        category_to_heading['bss'] = 'f_bss'
        category_to_heading['covalent'] = 'f_covalent'

        # load the .cif file
        structure = self.structure(pdb)

        # annotate pairwise interactions and extract loops
        interaction_to_triple_list, category_to_interactions, timerData, pair_to_data = annotate_nt_nt_in_structure(structure,categories)

        # actual categories found
        all_categories = list(set(categories.keys()) & set(category_to_interactions.keys()))

        pair_to_dictionary = {}

        for category in all_categories:
            if category == 'loops':
                continue

            heading = category_to_heading[category]

            for interaction in category_to_interactions[category]:

                if category == 'backbone':
                    if 'BPh' in interaction:
                        heading = 'f_bphs'
                    elif 'BR' in interaction:
                        heading = 'f_brbs'

                if category == 'basepair':
                    if "n" in interaction or "B" in interaction:
                        # no near basepairs, cWB, or cBW in basepair column
                        continue
                    # no lowercase to indicate geometry, like cWw or tsS
                    inter = interaction.replace("w","W").replace("s","S").replace("h","H")
                    # no basepair subcategories
                    inter = inter.replace("a","")
                else:
                    inter = interaction

                for triple in interaction_to_triple_list[interaction]:
                    u1 = triple[0]
                    u2 = triple[1]
                    pair = u1 + "+" + u2

                    if not pair in pair_to_dictionary:
                        pair_to_dictionary[pair] = {}

                    pair_to_dictionary[pair]['unit_id_1'] = u1
                    pair_to_dictionary[pair]['unit_id_2'] = u2
                    pair_to_dictionary[pair]['pdb_id'] = pdb
                    pair_to_dictionary[pair]['program'] = 'fr3d'
                    pair_to_dictionary[pair][heading] = inter

                    if 'f_crossing' in pair_to_dictionary[pair]:
                        if not pair_to_dictionary[pair]['f_crossing'] == triple[2]:
                            self.logger.info('Different crossing numbers for %s and %s interaction %s' % (u1,u2,interaction))
                    else:
                        pair_to_dictionary[pair]['f_crossing'] = triple[2]


        data_returned = False

        # return annotations
        for annotation in pair_to_dictionary.values():
            yield mod.UnitPairsInteractions2024(**annotation)
            data_returned = True

        # annotate self interactions and save in database
        self_interactions, error_message = annotate_self_base_backbone(structure)
        for d in self_interactions:
            if d['BPh']:
                unit_annotation_dict = {}
                unit_annotation_dict['pdb_id'] = pdb
                unit_annotation_dict['unit_id'] = d['unit_id']
                unit_annotation_dict['category'] = 'BPh'
                unit_annotation_dict['value'] = d['BPh']
                yield mod.UnitAnnotations(**unit_annotation_dict)
                data_returned = True
            if d['BR']:
                unit_annotation_dict = {}
                unit_annotation_dict['pdb_id'] = pdb
                unit_annotation_dict['unit_id'] = d['unit_id']
                unit_annotation_dict['category'] = 'BR'
                unit_annotation_dict['value'] = d['BR']
                yield mod.UnitAnnotations(**unit_annotation_dict)
                data_returned = True

        # save any loops that were extracted
        filename = os.path.join(DATA_FILE_DIRECTORY,'loops','%s_loops.txt' % pdb)
        with open(filename,'wt') as f:
            for full_loop in interaction_to_triple_list['loops']:
                a = full_loop['identifier']
                b = ",".join(full_loop['unit_ids'])
                c = ",".join(full_loop['border_indicators'])
                f.write('%s\t%s\t%s\n' % (a,b,c))

        if not data_returned:
            """
            # if we want to leave a note that we checked for interactions, we could try something like this:
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
            raise core.Skip("No new interactions to add for %s" % pdb)

