"""
Stage to search loops in one structure (or chains in the structure) against loops in other structures
(or chains in another structure) and write the matches to the loop_mapping table.
"""

import time

from cgi import print_arguments
import itertools as it
from operator import contains
import pymotifs.core as core
from pymotifs import models as mod
from pymotifs.utils import units
from pymotifs import utils
# from pymotifs.download import Downloader
# from pymotifs.pdbs.info import Loader as PdbLoader
# from pymotifs.units import Loader as UnitInfoLoder
from sqlalchemy import and_
from collections import defaultdict
from sqlalchemy import desc
from sqlalchemy import asc
from sqlalchemy import not_
from sqlalchemy import func
import annotations_for_chain_loops
from compare_and_cluster import list_all_chains_in
from external_motif_search import startup_external_dictionaries
from chain_to_rfam_family import read_equiv_class_csv_into_dict

from pymotifs.skip_files import SKIP
import sys

class Loader(core.Loader):
    merge_data = True
    mark = False
    allow_no_data = True
    # this stage is only run on existing chains, no dependencies needed
    dependencies = set([])

    def get_motif_atlas_loop_count(self):
        with self.session() as session:
            current_motif_atlas_release_query = session.query(mod.MlReleases.ml_release_id).\
                order_by(mod.MlReleases.date.desc()).\
                limit(1)
        for row in current_motif_atlas_release_query:
            current_motif_atlas_release = row.ml_release_id
        with self.session() as session:
            current_loop_count_in_motif_atlas = session.query(func.count(mod.MlLoops.loop_id)).\
                filter(mod.MlLoops.ml_release_id == current_motif_atlas_release).scalar()
        return current_loop_count_in_motif_atlas

    def to_process(self, pdbs, **kwargs):
        '''
        This method creates a list of things for this stage to work on.
        It returns a list of tuples.
        '''


        # where we store all of the chains to map
        # ex: list_to_work_on = [tuple([chain_to_annotate, representative_chains])]
        list_to_work_on = []
        chains_with_rfam_family = [] # useful for later, when doing chains w/o rfam families


        # grab all rfam families
        with self.session() as session:
            all_rfam_families_query = session.query(mod.ChainPropertyValue.value).\
                    filter(mod.ChainPropertyValue.property == "rfam_family").\
                    order_by(mod.ChainPropertyValue.value).\
                    distinct()

            all_rfam_families = []
            for row in all_rfam_families_query:
                all_rfam_families.append(row.value)


        # get pdbs that are xray diffraction and 3.5 angstroms or better
        with self.session() as session:
            trusted_pdb_query = session.query(mod.PdbInfo.pdb_id).\
                filter(mod.PdbInfo.experimental_technique.like("X-RAY%")).\
                filter(mod.PdbInfo.resolution <= 3.5)

            trusted_pdbs = []
            for row in trusted_pdb_query:
                trusted_pdbs.append(row.pdb_id)

        # this would be where to also append trusted cryo-em structures
        # trusted_cryo_em_query here
        # then append them to trusted_pdbs list

        # get from chain to equivalence class to representative
        chain_to_equiv_class_data = read_equiv_class_csv_into_dict(break_ifes = True)
        equiv_class_to_representatives = {}
        for chain, values in chain_to_equiv_class_data.items():
            if values['quality rank'] == 1:
                ec = values['equivalence class']

                if chain.split("|")[0] in trusted_pdbs:
                    if not equiv_class_to_representatives.get(ec):
                        equiv_class_to_representatives[ec] = [chain]
                    else:
                        if chain not in equiv_class_to_representatives[ec]:
                            equiv_class_to_representatives[ec] += [chain]


        # gather previous mappings, to avoid re-runs
        chain_to_previous_mappings = {} # will only have keys if that chain has been mapped to something before

        with self.session() as session:
            previous_mapping_query = session.query(mod.LoopMappingDone.unannotated_chain,
                    mod.LoopMappingDone.chain_mapped_to)

        for mapping in previous_mapping_query:
            if not chain_to_previous_mappings.get(mapping.unannotated_chain):
                chain_to_previous_mappings[mapping.unannotated_chain] = [mapping.chain_mapped_to]
            else:
                chain_to_previous_mappings[mapping.unannotated_chain].append(mapping.chain_mapped_to)


        # this is used as a trigger for when to run and when to skip previously run alignments
        motif_atlas_loop_count = self.get_motif_atlas_loop_count()

        # next line says "if motif atlas has 30 more loops than the last time many
        # chains were aligned", then run them all again
        if "RFAM CHAINS" in chain_to_previous_mappings: # safety first
            loop_counts = [int(i) for i in chain_to_previous_mappings["RFAM CHAINS"]]
            if motif_atlas_loop_count > max(loop_counts) + 30:
                rerun_aligned_chain_mappings = True
                self.logger.info("Rerunning aligned chain mappings due to 30 more loops in motif atlas")
            else:
                rerun_aligned_chain_mappings = False
        else:
            # this is to allow the first run of this pipeline stage to
            # build in chunks if it hits errors, instead of losing progress
            # DANGER: we NEED the final chunk of the first run to be at least
            #         100 chains long, else identifier row will need manually added
            #         to loop_mapping_done table. ie: ["RFAM CHAINS", current_motif_atlas_loop_count]
            rerun_aligned_chain_mappings = False

        # main loop
        for current_rfam in all_rfam_families:

            print("Looking for new chains in rfam family %s" % current_rfam)

            # get chains of given rfam family
            with self.session() as session:
                chains_of_rfam_query = session.query(mod.ChainPropertyValue.pdb_id,
                        mod.ChainPropertyValue.chain,
                        mod.ChainPropertyValue.property,
                        mod.ChainPropertyValue.value).\
                    filter(mod.ChainPropertyValue.property == "rfam_family").\
                    filter(mod.ChainPropertyValue.value == current_rfam)

                # find equiv classes from rfam family chain pool
                equiv_classes_in_this_rfam_family = []

                current_rfam_chains = []
                for row in chains_of_rfam_query:
                    # this line is useful later
                    chains_with_rfam_family.append("|".join([row.pdb_id, "1", row.chain]))
                    if row.pdb_id not in SKIP:
                        full_chain_id = "|".join([row.pdb_id, "1", row.chain])
                        current_rfam_chains.append(full_chain_id)
                        if chain_to_equiv_class_data.get(full_chain_id):
                            equiv_class_of_chain = chain_to_equiv_class_data[full_chain_id]["equivalence class"]
                            # print("chains: %s, ec: %s" % (full_chain_id, equiv_class_of_chain))
                            if equiv_class_of_chain not in equiv_classes_in_this_rfam_family:
                                # print(equiv_class_of_chain, equiv_class_to_representatives.get(equiv_class_of_chain), full_chain_id)
                                equiv_classes_in_this_rfam_family.append(equiv_class_of_chain)

            equiv_classes_in_this_rfam_family = sorted(equiv_classes_in_this_rfam_family)

            # equivalence to representatives
            # list of strings (strings are chains. IFEs are split and have the same quality rank for each chain)
            representatives_of_equiv_classes_in_this_rfam_family = []
            for equiv_class in equiv_classes_in_this_rfam_family:
                if equiv_class in equiv_class_to_representatives:
                    for chain in equiv_class_to_representatives[equiv_class]:
                        # with IFEs, there may be multiple chains here, include them all
                        if chain in current_rfam_chains:
                            representatives_of_equiv_classes_in_this_rfam_family.append(chain)

            # study RF02543
            # if current_rfam == "RF02543":
            #     #print("current_rfam_chains: %s" % current_rfam_chains)
            #     for crc in current_rfam_chains:
            #         if "5TBW" in crc:
            #             print("5TBW in current_rfam_chains: %s" % crc)
            #     print("equiv_class_to_representatives['NR_all_90506.65']: %s" % equiv_class_to_representatives["NR_all_90506.65"])
            #     print("equiv_classes_in_this_rfam_family: %s" % equiv_classes_in_this_rfam_family)
            #     print("representatives_of_equiv_classes_in_this_rfam_family: %s" % representatives_of_equiv_classes_in_this_rfam_family)

            #     time.sleep(15)

            arguments = [] # should get populated with tuples of (unannotated_chain, [trusted_chains])

            for unannotated_chain in current_rfam_chains:
                # if current_rfam == "RF02543":
                #     print("unannotated_chain: %s" % unannotated_chain)

                # insert logic here to detect 0 loop chains

                pdb_id, model_number, chain_id = unannotated_chain.split("|")

                if unannotated_chain not in representatives_of_equiv_classes_in_this_rfam_family:
                    for trusted_chain in representatives_of_equiv_classes_in_this_rfam_family:
                        if chain_to_previous_mappings.get(unannotated_chain):
                            if rerun_aligned_chain_mappings: # if we're forcing re-runs
                                arguments.append(tuple([unannotated_chain, representatives_of_equiv_classes_in_this_rfam_family]))
                                break
                            if trusted_chain not in chain_to_previous_mappings[unannotated_chain]:
                                arguments.append(tuple([unannotated_chain, representatives_of_equiv_classes_in_this_rfam_family]))
                                break # exits trusted_chain loop, to go to next "unannotated_chain in current_rfam_chains"
                        else: # if we have no mapping info about this chain yet, check if it has valid loops
                            with self.session() as session:
                                    # valid_loops_from_u_chain = session.query('''mod.loopInfo.loop_id, mod.loopInfo.loop_name''').\
                                    valid_loops_from_u_chain = session.query(mod.LoopInfo.loop_name).\
                                        filter(mod.LoopInfo.pdb_id == pdb_id).\
                                        filter(mod.LoopInfo.unit_ids.like("%" + unannotated_chain + "%")).count()

                            if valid_loops_from_u_chain > 0:
                                arguments.append((unannotated_chain, representatives_of_equiv_classes_in_this_rfam_family))
                            # else:
                                # print("skipped chain %s due to no valid loops" % unannotated_chain)
                                # sys.exit()
                            break # exits trusted_chain loop, to go to next "chain in current_rfam_chains"

            chains_added = set()
            for uc, reps in sorted(arguments):
                if not uc in chains_added:
                    chains_added.add(uc)
                    chains_added.update(reps)
                    list_to_work_on.append((uc,reps))
                    print("From Rfam family %s processing chains %s" % (current_rfam, uc))
                    self.logger.info("From Rfam family %s processing chains %s" % (current_rfam, uc))

            # list_to_work_on.extend(sorted(arguments))

        # insert a row in loop_mapping_done table as a way to keep track of how often to re-run everything
        # line indicates how many loops are in the motif atlas

        # instead, it seems that we could put the number as the second argument of the tuple, instead of a list of chains,
        # when there is no list of chains
        number_of_chains_with_rfam_families_to_be_run = len(list_to_work_on)

        if number_of_chains_with_rfam_families_to_be_run > 100:
            # the number 100 may need to be tuned. Im assuming no regular week will add more than 100
            # chains to rfam families.
            # but even if it does, we re-run early, and that's okay
            list_to_work_on.append(tuple(["RFAM CHAINS", motif_atlas_loop_count]))
            print("Planning to re-run Rfam chains")
            self.logger.info("Planning to re-run Rfam chains")
            # the above is saying "many alignments were run while there were X loops in the motif atlas"
            # if X is 30 more than the last time we re-ran mappings, then re-run them

        # gather chains without rfam families
        arguments = []
        with self.session() as session:
            all_chains_query = session.query(mod.ChainInfo.pdb_id,
                mod.ChainInfo.chain_name).\
                filter(mod.ChainInfo.entity_macromolecule_type.like("%polyribonucleotide%"))

            for row in all_chains_query:
                # create the chain id using model 1; later we will allow loops from
                # other models to match this pdb id and chain name
                full_chain_id = "|".join([row.pdb_id, "1", row.chain_name])
                if full_chain_id not in chains_with_rfam_family:
                    if full_chain_id in chain_to_previous_mappings:
                        # if the chain has been mapped to the motif atlas before,
                        # but now the motif atlas has at least 30 more loops than then,
                        # run again
                        # unannotated_chain: EXMPL|1|A, "MOTIF ATLAS"
                        # unannotated_chain: EXMPL|1|A, 3040
                        # ^ example of table rows if chain EXMPL|1|A were searched against
                        # 3040 loops in the motif atlas at that time
                        if "MOTIF ATLAS" in chain_to_previous_mappings[full_chain_id]:
                            highest_loop_count = 0
                            for loop_count in chain_to_previous_mappings[full_chain_id]:
                                if loop_count != "MOTIF ATLAS":
                                    if int(loop_count) > highest_loop_count:
                                        highest_loop_count = int(loop_count)
                            if motif_atlas_loop_count > highest_loop_count + 30:
                                arguments.append(tuple([full_chain_id, []]))
                                print("Mapping chain %s again; previously %s loops, now %d in Motif Atlas" % (full_chain_id,highest_loop_count,motif_atlas_loop_count))
                                self.logger.info("Mapping chain %s again; previously %s loops, now %d in Motif Atlas" % (full_chain_id,highest_loop_count,motif_atlas_loop_count))

                    else:
                        arguments.append(tuple([full_chain_id, []]))
                        print("Have not mapped chain %s before" % full_chain_id)
                        self.logger.info("Have not mapped chain %s before" % full_chain_id)

        list_to_work_on.extend(arguments)

        print("number of mappings to run: %s" % len(list_to_work_on))

        # sort in alphabetical order
        list_to_work_on = sorted(list_to_work_on, key=lambda x: x[0])

        # sort in reversed order
        list_to_work_on = sorted(list_to_work_on, key=lambda x: x[0], reverse=True)


        if len(list_to_work_on) == 0:
            self.logger.info("No loops to process")
            raise core.Skip("No loops to process")
        else:
            print(list_to_work_on[0:100])
            return list_to_work_on


    def has_data(self, pdb_dict, **kwargs):
        """
        This method can query the database after to_process but before data method
        to see if there is already data in the database, and if so, it returns True,
        and the data method does not have to work on this item.
        """

        return False

    def remove(self, pdb, **kwargs):

        """
        with self.session() as session:
            query = session.query(mod.LoopInfo).filter_by(pdb_id=pdb)
            ids = [result.loop_id for result in query]

        if not ids:
            return True

        with self.session() as session:
            return session.query(mod.LoopPositions).\
                filter(mod.LoopPositions.loop_id.in_(ids)).\
                delete(synchronize_session=False)
        """

        return True



    def data(self, input, **kwargs):
        """
        This method gets called on each item in the list returned by the
        to_process method, one after another.
        That will be one input at a time
        We will get one pdb identifier at a time, like pdb_id = '5HCQ'

        """

        # Adam start writing code here to deal with 5HCQ
        # Suggestion: use WinSCP to upload loops_and_strands.py to pymotifs/loops, don't commit
        # Suggestion: import * from loops_and_strands
        # Suggestion: hard code what you're searching for now
        # Add helper functions up above
        # Suggestion: do the searches


        # from annotations_for_chain_loops import main as mapper

        # print("I think chains_needing_annotated is: ", input[0])
        # print("I think suggested_annotated_chains are: ", input[1])
        chain_to_annotate = input[0]
        suggested_annotated_chains = input[1]

        motif_mapping_list = []
        formatted = {}

        print("#################################")
        print("working on annotating %s" % chain_to_annotate)

        # this is the core of the pipeline stage
        # motif_mapping_list is a list of dictionaries about successful matches
        # including loop id, query id, pdb ids, discrepancy

        max_size = 30  # largest loop size to consider, to save time by avoiding huge loops

        finished_loops = []
        with self.session() as session:
            # find loops in this PDB id that are already mapped
            # no need to search for them every time, since we seem to be getting multiple matches
            query = session.query(mod.LoopMapping.loop_id).\
                filter(mod.LoopMapping.pdb_id == chain_to_annotate.split("|")[0])

            for row in query:
                finished_loops.append(row.loop_id)

        motif_mapping_list = annotations_for_chain_loops.main(chain_to_annotate, suggested_annotated_chains, finished_loops, max_size)


        # import csv
        # with open("/usr/local/pipeline/hub-core/pymotifs/motif_atlas/5TBW-1-AR_vs_5TBW-1-1+7OSA-1-25S_then_vs_atlas.csv") as file:
        #     print("reading")
        #     reader = csv.DictReader(file)
        #     for row in reader:
        #         motif_mapping_list.append({'query_pdb_id': row['query_pdb_id'],
        #             'query_loop_id': row['query_loop_id'],
        #             'pdb_id': row['pdb_id'],
        #             'loop_id': row['loop_id'],
        #             'discrepancy': row['discrepancy'],
        #             'match_type' : row['match_type']})
                # print(line)
                # fields = line.split("\t")
                # formatted['query_pdb_id'] = fields[0]
                # formatted['query_loop_id'] = fields[1]
                # formatted['pdb_id'] = fields[2]
                # formatted['loop_id'] = fields[3]
                # formatted['discrepancy'] = fields[4]
                # formatted['match_type'] = fields[5]
                # motif_mapping_list.append(formatted)
                # motif_mapping_list.append({'query_pdb_id': fields[0],
                #     'query_loop_id': fields[1],
                #     'pdb_id': fields[2],
                #     'loop_id': fields[3],
                #     'discrepancy': fields[4]})
        if not suggested_annotated_chains:
            current_loop_count_in_motif_atlas = self.get_motif_atlas_loop_count()
            suggested_annotated_chains = ["MOTIF ATLAS", current_loop_count_in_motif_atlas]

        #self.logger.info("running on %s to %s" % (chain_to_annotate, suggested_annotated_chains))
        self.logger.info("length of motif_mapping_list: %d" % len(motif_mapping_list))


        if isinstance(suggested_annotated_chains, long):
            suggested_annotated_chains = [self.get_motif_atlas_loop_count()]
            # making it an iterable, for my indicator rows in loop_mapping_done

        # write what inputs have successfully run to the new database table loop_mapping_done
        for chain in suggested_annotated_chains:
            yield mod.LoopMappingDone(**{'unannotated_chain': chain_to_annotate, 'chain_mapped_to': chain})
            # print("skipping writing to LoopMappingDone")

        # yielding rows to the loop_mapping table if data is not empty
        for row in motif_mapping_list:
            # self.logger.info(row)
            # print(row) ###########################################################
            # row needs to be a dictionary about one motif mapping
            # keys are the column headers in the table
            # values are data values for the new row in the table

            # return data for the loop_mapping table
            yield mod.LoopMapping(**row)
        # else:
        #     # don't crash the pipeline, just move on with no data
        #     raise core.Skip("No matches for loops in %s found in %s" % (chain_to_annotate, suggested_annotated_chains))


