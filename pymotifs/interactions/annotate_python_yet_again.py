"""
Run Python code to annotate pairwise interactions
and store them in unit_pairs_annotations_2024, noting that program = 'fr3d'
2024 in the table name refers to the year when the parameters were set
Making a separate table will hopefully speed up queries for pairwise interactions
Next time there is a big change to the annotation program, make a new table
"""

import glob
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

    recompute = True  # force recompute of data; slow but thorough; deletes previous data

    # annotate structures with modified nts having no glycosidic center, because
    # these are the ones that may have been recently added to fr3d-python
    # This stage must be run before filling in the centers in order to work!
    annotate_modified = True
    annotate_modified = False

    modified_list = []
    modified_list = ['GUN','ADE','URA','TDR','A1I5K','CYT','A1CEB','A1IF9','A1L3O','A1BT7']

    process_all_by_rank = False  # process all files with a rank, by rank.  Bad idea.

    @property
    def table(self):
        return mod.UnitPairsInteractions2024

    def to_process(self, pdbs, **kwargs):
        """
        Return PDB ids which do not have python annotations in the table.
        Because of the nature of the annotations, even a file with only
        a few nucleotides should have at least one pairwise interaction.

        :param list pdbs: The list of pdb ids of interest on this run.
        :param dict kwargs: The keyword arguments which are ignored.
        :returns: A list of pdb ids to process.
        """

        if len(pdbs) == 1:
            return pdbs
        elif len(pdbs) < 99999:
            # the usual code to run
            # get (pdb_id,program) pairs in unit_pairs_interactions_2024 for these pdbs only; faster
            with self.session() as session:
                query = session.query(mod.UnitPairsInteractions2024.pdb_id,mod.UnitPairsInteractions2024.program).\
                filter(mod.UnitPairsInteractions2024.pdb_id.in_(pdbs)).\
                distinct()
                pdb_id_program_all = set([(result.pdb_id,result.program) for result in query])
            self.logger.info('Found %d pdb_id,program pairs' % len(pdb_id_program_all))
            pdbs_fr3d = set([pdb_id for pdb_id,program in pdb_id_program_all if program == 'fr3d'])
            need_to_annotate = sorted(pdbs_fr3d)
        elif self.annotate_modified:
            if len(self.modified_list) > 0:
                # process only files with the specified modified nucleotides
                with self.session() as session:
                    query = session.query(mod.UnitInfo.unit_id).\
                                    filter(mod.UnitInfo.unit_type_id.in_(['rna','dna'])).\
                                    filter(mod.UnitInfo.unit.in_(self.modified_list))
                    modified_na_units = set([r.unit_id for r in query])
                    self.logger.info('Found %d modified na units in list %s' % (len(modified_na_units),self.modified_list))
            else:
                # process all files with modified nucleotides, much slower
                with self.session() as session:
                    query = session.query(mod.UnitInfo.unit_id).\
                                    filter(mod.UnitInfo.unit_type_id.in_(['rna','dna'])).\
                                    filter(~mod.UnitInfo.unit.in_(['A','C','G','U','DA','DC','DG','DT','N','3DR']))
                    modified_na_units = set([r.unit_id for r in query])
                    self.logger.info('Found %d modified na units' % len(modified_na_units))

            # find all unit ids with glycosidic centers
            # takes almost 3 minutes
            with self.session() as session:
                query = session.query(mod.UnitCenters.unit_id).\
                               filter(mod.UnitCenters.name == 'glycosidic')
                glycosidic_units = set([r.unit_id for r in query])
                self.logger.info('Found %d units with glycosidic centers' % len(glycosidic_units))
            units_without_glycosidic = modified_na_units - glycosidic_units
            self.logger.info('Found %d modified RNA/DNA units without glycosidic')
            self.logger.info(str(sorted(units_without_glycosidic)))
            need_to_annotate = set([x.split("|")[0] for x in units_without_glycosidic])
        elif self.process_all_by_rank:
            # does not work as well as you might think, because some releases must
            # have given bad ranks to good structures
            # maybe order by nr_class_id, and then use the most recent rank
            # even then, obsolete files may have 0 as their most recent rank

            # maybe better yet, sort by rank 0, 1, 2, and by reverse nr_class_id

            # also misses many DNA structures!!!
            with self.session() as session:
                query = session.query(mod.NrClassRank.ife_id,mod.NrClassRank.rank,mod.NrClassRank.nr_class_id).\
                orderBy(mod.NrClassRank.nr_class_id).\
                distinct()

                pdb_to_rank = {}
                for result in query:
                    pdb = result.ife_id.split("|")[0]
                    rank = int(result.rank)
                    id = int(result.nr_class_id)
                    pdb_to_rank[pdb] = (rank,-id)

                pdb_rank = sorted(pdb_to_rank.items(), key=lambda x : (x[1][0],x[1][1],x[0]))

                self.logger.info(pdb_rank)

                need_to_annotate = [x[0] for x in pdb_rank]
        else:
            # get all (pdb_id,program) pairs in unit_pairs_interactions_2024 table
            with self.session() as session:
                query = session.query(mod.UnitPairsInteractions2024.pdb_id,mod.UnitPairsInteractions2024.program).\
                distinct()
                pdb_id_program_all = set([(result.pdb_id,result.program) for result in query])
            self.logger.info('Found %d pdb_id,program pairs' % len(pdb_id_program_all))
            pdbs_fr3d = set([pdb_id for pdb_id,program in pdb_id_program_all if program == 'fr3d'])
            need_to_annotate = sorted(pdbs_fr3d)


        if False:
            # temporary to get all pdb files in unit_info
            # (unit ids will exist and not cause foreign key problems)
            with self.session() as session:
                query = session.query(mod.UnitInfo.pdb_id).\
                distinct()
                pdbs = [result.pdb_id for result in query]
            print('Found %d pdb files in unit_info' % len(pdbs))

            # get a list of .cif files in /usr/local/pipeline/hub-core/
            cif_files = os.listdir('/usr/local/pipeline/hub-core/cif_files/')
            cif_files = [f for f in cif_files if f.endswith('.cif.gz')]
            cif_files = [f.replace('.cif.gz','') for f in cif_files]

            pdbs_fr3d = sorted(set(pdbs) & set(cif_files))

        if False:
            # get all pdb files that are the basis for the motif atlas
            # won't take as long to process those
            nr_release_id = '3.352'  # for release 3.88
            nr_release_id = '3.360'  # for release 3.90
            class_start = 'NR'
            from pymotifs.constants import MOTIF_ALLOWED_METHODS
            from pymotifs.constants import MOTIF_RESOLUTION_CUTOFF

            motif_atlas_pdbs = set()
            with self.session() as session:
                classranks = mod.NrClassRank
                classes = mod.NrClasses
                ifes = mod.IfeInfo
                pdbs = mod.PdbInfo
                query = session.query(classranks,pdbs.pdb_id).\
                    join(classes, classes.name == classranks.nr_class_name).\
                    join(ifes, ifes.ife_id == classranks.ife_id).\
                    join(pdbs, pdbs.pdb_id == ifes.pdb_id).\
                    filter(classranks.rank == 0).\
                    filter(classes.nr_release_id == nr_release_id).\
                    filter(classes.resolution == MOTIF_RESOLUTION_CUTOFF).\
                    filter(classes.name.like('%s%%' % class_start)).\
                    filter(ifes.new_style == True).\
                    filter(pdbs.experimental_technique.in_(MOTIF_ALLOWED_METHODS))

                for row in query:
                    motif_atlas_pdbs.add(row.pdb_id)

                print('Found %d motif atlas files' % len(motif_atlas_pdbs))

        # # H.m. structures to practice on
        # need_to_annotate = ['1FFK','1S72','4V9F']

        # # annotate all of the files; do every other file starting at 0
        # need_to_annotate = need_to_annotate[0::2]
        # # annotate all of the files; do every other file starting at 1
        # need_to_annotate = need_to_annotate[1::2]

        # # annotate all of the structures needed for the next motif atlas release
        # need_to_annotate = motif_atlas_pdbs

        # # if we already did the motif atlas files, just do the rest
        # need_to_annotate = sorted(set(pdbs_fr3d) - motif_atlas_pdbs)
        # self.logger.info('Processing everything other than 3.360 files for motif atlas')

        print('Found %d files to annotate' % len(need_to_annotate))

        if len(need_to_annotate) == 0:
            raise core.Skip("No PDB files need pairwise interactions annotated with fr3d-python")
        else:
            return need_to_annotate


    def query(self, session, pdb):
        return session.query(mod.UnitPairsInteractions2024).\
                filter_by(pdb_id=pdb).\
                filter_by(program='fr3d')


    def convert_symmetry(self,u,old_prefix,new_prefix):
        # convert ASM_ to P_ for saving in the database
        if old_prefix:
            f = u.split("|")
            if len(f) == 9:
                if f[8].startswith(old_prefix):
                    f[8] = f[8].replace(old_prefix,new_prefix)
                    u = "|".join(f)
        return u


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

        # naming convention for symmetry operators in pairwise annotations
        annotation_symmetry = set()
        for interaction, triple_list in interaction_to_triple_list.items():
            if interaction == 'loops':
                continue
            for u1, u2, crossing in triple_list:
                f1 = u1.split("|")
                if len(f1) == 9:
                    if f1[8].startswith('ASM'):
                        annotation_symmetry.add('ASM')
                    elif f1[8].startswith('P'):
                        annotation_symmetry.add('P')
                f2 = u2.split("|")
                if len(f2) == 9:
                    if f2[8].startswith('ASM'):
                        annotation_symmetry.add('ASM')
                    elif f2[8].startswith('P'):
                        annotation_symmetry.add('P')

        self.logger.info('annotation_symmetry is %s' % annotation_symmetry)

        # naming convention for symmetry operators in existing entry like 6GV4
        old_prefix = ""
        new_prefix = ""
        if 'ASM' in annotation_symmetry:
            # check to see if units in the database have P_ symmetry instead
            with self.session() as session:
                query = session.query(mod.UnitInfo.sym_op).\
                filter(mod.UnitInfo.pdb_id == pdb).\
                distinct()
                database_symmetry = set([result.sym_op.split("_")[0] for result in query])

            self.logger.info('database_symmetry is %s' % database_symmetry)

            if 'P' in database_symmetry and not 'ASM' in database_symmetry:
                old_prefix = "ASM_"
                new_prefix = "P_"

                self.logger.info('Changing unit ids in new annotations from %s to %s' % (old_prefix,new_prefix))

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

                    if new_prefix:
                        u1 = self.convert_symmetry(u1,old_prefix,new_prefix)
                        u2 = self.convert_symmetry(u2,old_prefix,new_prefix)

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
                unit_annotation_dict['unit_id'] = self.convert_symmetry(d['unit_id'],old_prefix,new_prefix)
                unit_annotation_dict['category'] = 'BPh'
                unit_annotation_dict['value'] = d['BPh']
                yield mod.UnitAnnotations(**unit_annotation_dict)
                data_returned = True
            if d['BR']:
                unit_annotation_dict = {}
                unit_annotation_dict['pdb_id'] = pdb
                unit_annotation_dict['unit_id'] = self.convert_symmetry(d['unit_id'],old_prefix,new_prefix)
                unit_annotation_dict['category'] = 'BR'
                unit_annotation_dict['value'] = d['BR']
                yield mod.UnitAnnotations(**unit_annotation_dict)
                data_returned = True

        # save any loops that were extracted
        filename = os.path.join(DATA_FILE_DIRECTORY,'loops','%s_loops.txt' % pdb)
        with open(filename,'wt') as f:
            for full_loop in interaction_to_triple_list['loops']:
                a = full_loop['identifier']
                if old_prefix:
                    new_ids = [self.convert_symmetry(u,old_prefix,new_prefix) for u in full_loop['unit_ids']]
                    b = ",".join(new_ids)
                else:
                    b = ",".join(full_loop['unit_ids'])
                c = ",".join(full_loop['border_indicators'])
                f.write('%s\t%s\t%s\n' % (a,b,c))

        if data_returned:
            # delete pickle files with pairwise interactions
            # so they will be exported again the next time
            for f in glob.glob("/usr/local/pipeline/hub-core/data/pairs/%s_NA_pairs.pickle" % pdb):
                os.remove(f)
        else:
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

