"""
Stage to search loops in one structure (or chains in the structure) against loops in other structures
(or chains in another structure) and write the matches to the loop_mapping table.
"""

from sqlalchemy import func

from pymotifs import core
from pymotifs import models as mod

# from compare_and_cluster import list_all_chains_in
# from external_motif_search import startup_external_dictionaries
from pymotifs.motif_atlas.chain_to_rfam_family import read_equiv_class_csv_into_dict
from pymotifs.motif_atlas.annotations_for_chain_loops import main as annotation_main
from pymotifs.motif_atlas.annotations_for_chain_loops import update_motif_atlas_files

from pymotifs.skip_files import SKIP

class Loader(core.SimpleLoader):
    merge_data = True
    mark = False
    allow_no_data = True
    # this stage is only run on existing chains, no dependencies needed
    dependencies = set([])

    def get_motif_atlas_loop_count(self,current_motif_atlas_release):
        # with self.session() as session:
        #     current_motif_atlas_release_query = session.query(mod.MlReleases.ml_release_id).\
        #         order_by(mod.MlReleases.date.desc()).\
        #         limit(1)
        # for row in current_motif_atlas_release_query:
        #     current_motif_atlas_release = row.ml_release_id

        with self.session() as session:
            current_loop_count_in_motif_atlas = session.query(func.count(mod.MlLoops.loop_id)).\
                filter(mod.MlLoops.ml_release_id == current_motif_atlas_release).scalar()
        return current_loop_count_in_motif_atlas


    def to_process(self, pdbs, **kwargs):
        """
        This method creates a list of things for this stage to work on.
        It returns a list of tuples.
        """

        # control what representative set release we compare to
        # every week, we should probably compare to the previous week
        # but if we miss something, it's nice to be able to go back and compare to earlier ones
        release = 'current'
        release = '3.362'
        release = 'previous'   # use each week, to compare to last week's representatives

        # where we store all of the chains to map
        # ex: list_to_work_on = [tuple([chain_to_annotate, representative_chains])]
        list_to_work_on = []
        chains_with_rfam_family = [] # useful for later, when doing chains w/o rfam families

        # get all rfam families that have PDB chains mapped to them
        all_rfam_families = []
        with self.session() as session:
            all_rfam_families_query = session.query(mod.ChainPropertyValue.value).\
                    filter(mod.ChainPropertyValue.property == "rfam_family").\
                    order_by(mod.ChainPropertyValue.value).\
                    distinct()

            for row in all_rfam_families_query:
                all_rfam_families.append(row.value)

        self.logger.info('Found %d Rfam families with PDB chains' % len(all_rfam_families))

        # get pdbs that are x-ray diffraction and 3.5 angstroms or better
        # these are the ones that the motif atlas is based on
        trusted_pdbs = []
        with self.session() as session:
            trusted_pdb_query = session.query(mod.PdbInfo.pdb_id).\
                filter(mod.PdbInfo.experimental_technique.like("X-RAY%")).\
                filter(mod.PdbInfo.resolution <= 3.5)

            for row in trusted_pdb_query:
                trusted_pdbs.append(row.pdb_id)

        self.logger.info('Found %d trusted PDBs' % len(trusted_pdbs))

        # this would be where to also append trusted cryo-em structures
        # then append them to trusted_pdbs list

        # alternatively, process the list of loops in the release of the motif atlas
        # and use those PDBs

        # map from chain/ife to equivalence class to representative
        chain_to_equiv_class_data = read_equiv_class_csv_into_dict(release = release, break_ifes = True)
        equiv_class_to_representatives = {}
        for chain, values in chain_to_equiv_class_data.items():
            # keep chains that are the representative of their set
            if values['quality rank'] == 1:
                ec = values['equivalence class']

                if chain.split("|")[0] in trusted_pdbs:
                    if not equiv_class_to_representatives.get(ec):
                        equiv_class_to_representatives[ec] = [chain]
                    else:
                        if chain not in equiv_class_to_representatives[ec]:
                            equiv_class_to_representatives[ec] += [chain]

        # gather previous mappings, to avoid re-runs
        # key will be a chain, value will be a list of chains it was mapped to
        # will only have keys if that chain has been mapped to something before
        chain_to_previous_mappings = {}

        chain_to_motif_atlas_size = {}
        with self.session() as session:
            previous_mapping_query = session.query(mod.LoopMappingDone.unannotated_chain,
                    mod.LoopMappingDone.chain_mapped_to)

            for mapping in previous_mapping_query:
                unannotated_chain = mapping.unannotated_chain
                chain_mapped_to = mapping.chain_mapped_to

                if 'ATLAS ' in chain_mapped_to:
                    motif_atlas_size = int(chain_mapped_to.split(" ")[1])
                    if unannotated_chain in chain_to_motif_atlas_size:
                        # keep track of the largest number of motif atlas loops compared to
                        chain_to_motif_atlas_size[unannotated_chain] = max(motif_atlas_size,chain_to_motif_atlas_size[unannotated_chain])
                    else:
                        # first mention of how many loops were in the motif atlas
                        chain_to_motif_atlas_size[unannotated_chain] = motif_atlas_size
                else:
                    if chain_to_previous_mappings.get(unannotated_chain):
                        chain_to_previous_mappings[unannotated_chain].append(mapping.chain_mapped_to)
                    else:
                        chain_to_previous_mappings[unannotated_chain] = [mapping.chain_mapped_to]

        # find the most recent motif atlas release, create queries and search spaces,
        # save as .pickle file so it's quick to load for each new chain
        release = update_motif_atlas_files()

        # how many loops are in the current motif atlas
        motif_atlas_loop_count = self.get_motif_atlas_loop_count(release)
        self.logger.info('There are %d loops in the current motif atlas' % motif_atlas_loop_count)

        # loop over Rfam families, find the chains mapped to each Rfam family
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

            # equivalence class to representatives
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
                        if len(pdbs) == 1 and pdbs[0] in unannotated_chain:
                            # we are specifically running on one pdb id, so process it
                            arguments.append((unannotated_chain, representatives_of_equiv_classes_in_this_rfam_family))
                            break
                        elif chain_to_previous_mappings.get(unannotated_chain):
                            # we have already annotated this chain before, but maybe we should try it again
                            if motif_atlas_loop_count > chain_to_motif_atlas_size.get(unannotated_chain,0) + 50:
                                arguments.append((unannotated_chain, representatives_of_equiv_classes_in_this_rfam_family))
                                self.logger.info("Mapping chain %s again; previously %s loops, now %d in Motif Atlas" % (unannotated_chain,chain_to_motif_atlas_size.get(unannotated_chain,0),motif_atlas_loop_count))
                                break
                            if trusted_chain not in chain_to_previous_mappings[unannotated_chain]:
                                arguments.append((unannotated_chain, representatives_of_equiv_classes_in_this_rfam_family))
                                break # exits trusted_chain loop, to go to next "unannotated_chain in current_rfam_chains"
                        else:
                            # if we have no mapping info about this chain yet, check if it has valid loops
                            with self.session() as session:
                                valid_loops_from_u_chain = session.query(mod.LoopInfo.loop_name).\
                                    join(mod.LoopPositions, mod.LoopInfo.loop_id == mod.LoopPositions.loop_id).\
                                    filter(mod.LoopInfo.pdb_id == pdb_id).\
                                    filter(mod.LoopPositions.unit_id.like(unannotated_chain + "%")).count()

                            if valid_loops_from_u_chain > 0:
                                self.logger.info("Mapping chain %s because it is not in chain_to_previous_mappings" % (unannotated_chain))
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
                    list_to_work_on.append((uc,reps,release,motif_atlas_loop_count))
                    print("From Rfam family %s processing chains %s" % (current_rfam, uc))
                    self.logger.info("From Rfam family %s processing chains %s" % (current_rfam, uc))

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
                    highest_loop_count = chain_to_motif_atlas_size.get(full_chain_id,0)
                    if len(pdbs) == 1 and pdbs[0] in full_chain_id:
                        # we are specifically running on one pdb id, so process it
                        arguments.append((full_chain_id, [], release, motif_atlas_loop_count))
                        print("Mapping chain %s; previously %s loops, now %d in Motif Atlas" % (full_chain_id,highest_loop_count,motif_atlas_loop_count))
                        self.logger.info("Mapping chain %s; previously %s loops, now %d in Motif Atlas" % (full_chain_id,highest_loop_count,motif_atlas_loop_count))
                    elif full_chain_id in chain_to_previous_mappings or full_chain_id in chain_to_motif_atlas_size:
                        # if the chain has been mapped to the motif atlas before,
                        # but now the motif atlas has at least 30 more loops than then,
                        # run again

                        # we have already annotated this chain before, but maybe we should try it again
                        if motif_atlas_loop_count > highest_loop_count + 50:
                            arguments.append((full_chain_id, [], release, motif_atlas_loop_count))
                            print("Mapping chain %s again; previously %s loops, now %d in Motif Atlas" % (full_chain_id,highest_loop_count,motif_atlas_loop_count))
                            self.logger.info("Mapping chain %s again; previously %s loops, now %d in Motif Atlas" % (full_chain_id,highest_loop_count,motif_atlas_loop_count))

                    else:
                        arguments.append((full_chain_id, [], release, motif_atlas_loop_count))
                        print("Have not mapped chain %s before" % full_chain_id)
                        self.logger.info("Have not mapped chain %s before" % full_chain_id)

        list_to_work_on.extend(arguments)

        new_list = []
        for entry in list_to_work_on:
            # remove problematic entries
            # 4V3P has thousands of loops!
            if '4V3P' in entry[0]:
                continue

            if len(pdbs) == 1 and not pdbs[0] in entry[0]:
                continue

            new_list.append(entry)

        list_to_work_on = new_list

        print("number of mappings to run: %s" % len(list_to_work_on))

        # sort chains in alphabetical order
        list_to_work_on = sorted(list_to_work_on, key=lambda x: x[0])

        # sort chains in reverse alphabetical order, to process the most recent PDB files first
        list_to_work_on = sorted(list_to_work_on, key=lambda x: x[0], reverse=True)


        """
        # target very specific chains and mappings

        # 8C3A released 2024-01-10, just after release 3.79
        release = update_motif_atlas_files('3.79')
        list_to_work_on = [('8C3A|1|1',['5TBW|1|1'],'3.79')]
        list_to_work_on = [('8C3A|1|CM',['4V88|1|A6'],'3.79')]

        # 8VTW released 2024-08-07, just after release 3.86
        release = update_motif_atlas_files('3.86')
        list_to_work_on = [('8VTW|1|1A',['5J7L|1|DA','7RQB|1|1A','7A0S|1|X','4WF9|1|X'],'3.86')]
        """

        if len(list_to_work_on) == 0:
            raise core.Skip("No loops to process")
        else:
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
        input is a tuple like ('9DFD|1|1a', ['5J7L|1|AA', '4LFB|1|A', '6CZR|1|1a'], 3.88, 4032)
        """

        # if input[0] == 'RFAM CHAINS':
        #     # special case, to record the size of the motif atlas
        #     yield mod.LoopMappingDone(**{'unannotated_chain': input[0], 'chain_mapped_to': input[1]})
        #     return

        chain_to_annotate = input[0]
        pdb_id = chain_to_annotate.split("|")[0]
        suggested_annotated_chains = input[1]           # could be empty, when no Rfam family
        release = input[2]
        current_loop_count_in_motif_atlas = input[3]

        # motif_mapping_list is a list of dictionaries about successful matches
        # including loop id, query id, pdb ids, discrepancy
        motif_mapping_list = []

        max_size = 30  # largest loop size to consider, to save time by avoiding huge loops

        finished_loops = []
        with self.session() as session:
            # find loops in this PDB id that are already mapped
            query = session.query(mod.LoopMapping.loop_id).\
                filter(mod.LoopMapping.pdb_id == pdb_id)

            for row in query:
                finished_loops.append(row.loop_id)

        self.logger.info('Starting to annotate loops in %s' % chain_to_annotate)
        print('')
        print('Starting to annotate loops in %s' % chain_to_annotate)

        motif_mapping_list = annotation_main(chain_to_annotate, suggested_annotated_chains, finished_loops, max_size, release)

        #self.logger.info("running on %s to %s" % (chain_to_annotate, suggested_annotated_chains))
        self.logger.info("length of motif_mapping_list: %d" % len(motif_mapping_list))

        # write what inputs have successfully run to the new database table loop_mapping_done
        for chain in suggested_annotated_chains:
            yield mod.LoopMappingDone(**{'unannotated_chain': chain_to_annotate, 'chain_mapped_to': chain})
            # print("skipping writing to LoopMappingDone")

        # note how many motif atlas loops this chain has (potentially) been compared against
        yield mod.LoopMappingDone(**{'unannotated_chain': chain_to_annotate, 'chain_mapped_to': "ATLAS %d" % current_loop_count_in_motif_atlas})

        # yield rows to the loop_mapping table if data is not empty
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
