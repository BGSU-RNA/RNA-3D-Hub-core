'''
the main purpose of this script is to find annotations for
loops, given a chain
by Adam Coger, acoger@bgsu.edu 04/24/2023
'''

import sys
import os
import time
sys.path.append('search')

from fr3d_configuration import DATAPATHUNITS, DATAPATHRESULTS, DATAPATHATLAS, DATAPATHLOOPS, SERVER
from external_motif_search import make_list_of_all_loops_in, startup_external_dictionaries
from compare_and_cluster import load_motif_atlas_release, create_queries_and_search_spaces, save_motif_atlas_already_setup, load_motif_atlas_already_setup, run_aligned_searches_first, all_against_all_searches, load_previous_search_results_one_pdb
import csv # temporary
from compare_motifs_of_two_structures import write_output_csv


def get_list_of_relevant_motif_atlas_loops(types = ["HL", "IL"], release = "current"):
    loop_ids = []
    motif_list = []

    # loop_type refers to name from number of strands
    for loop_type in types:
        print()
        print("trying to load: %s" % loop_type)
        print()
        motif_list.extend(load_motif_atlas_release(loop_type = loop_type, release = release, path = DATAPATHATLAS))

    for motif in motif_list:
        loop_ids.extend(motif['alignment'].keys())

    return(loop_ids)


def main(chain_id = "6ZMI|1|S2", chains_for_comparing_against = [], finished_loops = [], max_size = None):
    # chain_id is allowed to be an IFE

    # if someone input a string as the second variable
    # this prevents indexing going to characters
    if isinstance(chains_for_comparing_against, str):
        chains_for_comparing_against = [chains_for_comparing_against]
    if isinstance(chains_for_comparing_against, long):
        return()


    # preparing chain loops for FR3D / AvA searches
    loop_ids_of_chain = make_list_of_all_loops_in(structure = chain_id, max_size = max_size)

    loop_ids_of_chain = sorted(set(loop_ids_of_chain) - set(finished_loops))

    print("Loops in chain %s are %s" % (chain_id,sorted(loop_ids_of_chain)))

    loops_of_chain = startup_external_dictionaries(loop_ids_of_chain)
    chain_as_queries, chain_as_search_spaces, chain_flanking_bp_queries = create_queries_and_search_spaces(loops_of_chain)

    if len(chain_as_search_spaces) == 0:
        # skip a lot of work if chain has issues
        return()

    # I should find a way to detect missing nts...
    problems = ["HL_6ZMI_059", "J3_6ZMI_003"]
    for problem in problems:
        if(problem in chain_as_queries):
            del chain_as_queries[problem]
            del chain_as_search_spaces[problem]
            del chain_flanking_bp_queries[problem]

    release = "current"
    # gather all motif atlas loops
    files_already_exist = bool(os.path.isfile(os.path.join(DATAPATHATLAS, "%s_atlas_queries.pickle" % release)) and
        os.path.isfile(os.path.join(DATAPATHATLAS, "%s_atlas_search_spaces.pickle" % release)) and
        os.path.isfile(os.path.join(DATAPATHATLAS, "%s_atlas_flanking_bp_queries.pickle" % release)))
    one_month = 30*24*60*60 # in seconds
    if files_already_exist:
        # see if their file modification date + 1 month is greater than today
        files_are_up_to_date = bool(os.path.getmtime(os.path.join(DATAPATHATLAS, "%s_atlas_queries.pickle" % release)) + one_month > time.time() and
            os.path.getmtime(os.path.join(DATAPATHATLAS, "%s_atlas_search_spaces.pickle" % release)) + one_month > time.time() and
            os.path.getmtime(os.path.join(DATAPATHATLAS, "%s_atlas_flanking_bp_queries.pickle" % release)) + one_month > time.time())
    else:
        files_are_up_to_date = False # just preventing an error
    if not files_already_exist or not files_are_up_to_date:
        loop_types = ["HL", "IL"] # , "J3"] # j3 leads to a php error
        delete_me = ""
        motif_atlas_loop_ids = get_list_of_relevant_motif_atlas_loops(types = loop_types)
        loops_of_atlas = startup_external_dictionaries(motif_atlas_loop_ids)
        atlas_queries, atlas_search_spaces, atlas_flanking_bp_queries = create_queries_and_search_spaces(loops_of_atlas)
        save_motif_atlas_already_setup(release, atlas_queries, atlas_search_spaces, atlas_flanking_bp_queries)
        print("updated motif atlas files")
    else:
        atlas_queries, atlas_search_spaces, atlas_flanking_bp_queries = load_motif_atlas_already_setup(release)
        print("loaded motif atlas files")
    # loop over chains to compare against
    # instantiate accumulator dicts
    double_filtered_results = {}
    aligned_results = {}
    unmatched_search_spaces = {} # this is the overlap of leftover_search_spaces from ALL annotated chains for comparing against
    search_results = {} # populated from motif atlas geometric searches
    for annotated_chain_for_comparing in chains_for_comparing_against:
        # print("in comparison loop. checking %s against %s" % (chain_id, annotated_chain_for_comparing))
        # gather loops from this annotated chain
        loop_ids_of_annotated_chain = make_list_of_all_loops_in(structure = annotated_chain_for_comparing)
        loops_of_annotated_chain = startup_external_dictionaries(loop_ids_of_annotated_chain)
        annotated_chain_as_queries, annotated_chain_as_search_spaces, annotated_chain_flanking_bp_queries = create_queries_and_search_spaces(loops_of_annotated_chain)


        # look for alignment and run 1 to 1 FR3D searches on aligned loops
        priority_results, leftover_search_spaces = run_aligned_searches_first(queries = annotated_chain_as_queries,
            search_spaces = chain_as_search_spaces, flanking_bp_queries = annotated_chain_flanking_bp_queries,
            force_alignment = True)

        aligned_results.update(priority_results)
        if not unmatched_search_spaces:
            unmatched_search_spaces = leftover_search_spaces
        else:
            # find intersection/overlap of search spaces that had no matches
            if sys.version_info[0] < 3: # running python 2
                loops_in_common = unmatched_search_spaces.viewkeys() & leftover_search_spaces.viewkeys()
            else: # python 3+
                loops_in_common = unmatched_search_spaces.keys() & leftover_search_spaces.keys()
            unmatched_search_spaces = {loop: chain_as_search_spaces[loop] for loop in loops_in_common}

    if not unmatched_search_spaces:
        unmatched_search_spaces = chain_as_search_spaces

    # check if we eliminated any loop types to shorten the next step
    loop_types = set()
    # if unmatched_search_spaces:
    #     for loop_id in unmatched_search_spaces.keys():
    #         loop_types.add(loop_id.split("_")[0])
    # else:
    for loop_id in chain_as_search_spaces:
        loop_types.add(loop_id.split("_")[0])

    # run AvA searches from motif atlas queries to unmatched search spaces
    for prefix in loop_types:
        loop_type_queries = {loop: query_info for loop, query_info in atlas_queries.items() if loop.startswith(prefix)}
        loop_type_search_spaces = {loop: query_info for loop, query_info in unmatched_search_spaces.items() if loop.startswith(prefix)}
        print('searching for %s matches from %d loops' % (prefix, len(loop_type_queries)))

        if len(loop_type_search_spaces) > 0:
            loop_type_results = all_against_all_searches(loop_type = prefix,
                queries = loop_type_queries, search_spaces = loop_type_search_spaces,
                flanking_bp_queries = atlas_flanking_bp_queries, load_saved_searches = True,
                save_path = os.path.join(DATAPATHRESULTS, "for_annotation"))
            search_results.update(loop_type_results)
        if len(loop_type_search_spaces) < 2:
            # if there are no search spaces fed to all_against_all, it will not load prev_results
                print("in force load for %s loops" % prefix)
                results = load_previous_search_results_one_pdb(loop_type = prefix, pdb_id = chain_id.split("|")[0], path = os.path.join(DATAPATHRESULTS, "for_annotation"))
                for pair, info in results:
                    if pair[0] in atlas_queries and pair[1] in chain_as_search_spaces:
                        if not isinstance(results[(query_id, search_space_id)],int): # says "if not disqualified"
                            # print("added a match!")
                            search_results[(query_id,search_space_id)] = results[(query_id,search_space_id)]
        # print(" %s : %d" % (prefix, len(loop_type_search_spaces)))

    # cutting out disqualified and high discrepancy results
    filtered_results = {}
    exhaustive_discrepancy_cutoff = .4
    for pair, info in search_results.items():
        if pair[1] in chain_as_search_spaces:
            if info[0]['dq'] == []:
                if info[0]['discrepancy'] <= exhaustive_discrepancy_cutoff:
                    filtered_results[pair] = info
    print("unannotated chain: %s contains %d loops" % (chain_id, len(chain_as_search_spaces)))
    # print("there were %d homologous matches, before filtering" % len(aligned_results))
    # print("there were %d geometric matches, before filtering" % len(filtered_results))
    filtered_results.update(aligned_results) # to dodge exhaustive cutoff above, we add homologous matches here

    # this shows that we end up only with BEST MATCHES and don't prioritize homologous?!
    # for search_space in chain_as_search_spaces:
    #     best_discrepancy = 100
    #     for pair, info in filtered_results.items():
    #         if search_space in pair:
    #             if info[0]['discrepancy'] < best_discrepancy:
    #                 best_discrepancy = info[0]['discrepancy']
    #                 best_pair = pair
    #                 best_info = info
    #     if best_discrepancy != 100:
    #         double_filtered_results[best_pair] = best_info

    # edited to prioritize homologous
    for search_space in chain_as_search_spaces:
        best_geo_discrepancy = 100
        best_homo_discrepancy = 100
        for pair, info in filtered_results.items():
            if search_space in pair:
                if info[0]['match_type'] == "geometric":
                    if info[0]['discrepancy'] < best_geo_discrepancy:
                        best_geo_discrepancy = info[0]['discrepancy']
                        best_geo_pair = pair
                        best_geo_info = info
                if info[0]['match_type'] == "homologous":
                    if info[0]['discrepancy'] < best_homo_discrepancy:
                        best_homo_discrepancy = info[0]['discrepancy']
                        best_homo_pair = pair
                        best_homo_info = info
        if best_geo_discrepancy != 100 or best_homo_discrepancy != 100:
            if best_homo_discrepancy != 100:
                double_filtered_results[best_homo_pair] = best_homo_info
                # print("best homologous match: q %s, ss %s" % (best_homo_pair[0], best_homo_pair[1]))
            elif best_geo_discrepancy != 100: # only record geometric match if we lack a homologous one
                double_filtered_results[best_geo_pair] = best_geo_info
                # print("best geometric match: q %s, ss %s" % (best_geo_pair[0], best_geo_pair[1]))


    # reprocess results to fit LOOP_MAPPING table on server
    # fields: query_pdb_id, query_loop_id, pdb_id, loop_id, discrepancy, match_type
    results_in_list_of_dicts_format = []
    for searched_pair, details in double_filtered_results.items():
        # motif_atlas_size = len(motif_atlas_queries)
        q_pdb_id = details[0]["query_id"].split("_")[1]
        pdb_id = details[0]["search_space_id"].split("_")[1]
        results_in_list_of_dicts_format.append({"query_pdb_id": q_pdb_id, "query_loop_id": details[0]["query_id"],
            "pdb_id": pdb_id, "loop_id": details[0]["search_space_id"], "discrepancy": details[0]["discrepancy"],
            "match_type": details[0]["match_type"]})
    print("we have %d rows to write to loop_mapping table" % len(results_in_list_of_dicts_format))
    # writing a csv file in a different format for html comparison
    # write_output_csv(results = double_filtered_results, filename = "%s_vs_%s_then_vs_atlas.csv" % (chain_id, annotated_chain_for_comparing))
    # if sys.version_info[0] >= 3:
    #     write_output_csv(results = double_filtered_results, filename = "%s_vs_%s_then_vs_atlas.csv" % (chain_id, "+".join(chains_for_comparing_against)))

    return(results_in_list_of_dicts_format)



if __name__ == "__main__":

    # unannotated_chain = '6ZMI|1|L5'
    # annotated_chains_to_align = ['5TBW|1|1']

    # unannotated_chain = '5TBW|1|AR'
    # annotated_chains_to_align = ['5TBW|1|1', '7OSA|1|25S']

    # testing for the awful IL_4WFB_107 to IL_4WF9_111 FR3D Search
    unannotated_chain = '4WFB|1|Y'
    annotated_chains_to_align = ['4WF9|1|Y']
    # annotated_chains_to_align = []

    results = main(unannotated_chain, annotated_chains_to_align)
    print("\a")
    print("\a")
    print("\a")

    # file writing moved inside of main()
    # filename = "%s_vs_%s_then_vs_atlas.csv" % (unannotated_chain.replace("|", "-"),
    #     annotated_chain_to_align.replace("|", "-"))

    # for the below line to work, main() needs to
    # return "double_filtered_results" instead of the
    # current format used to fill SQL table
    # write_output_csv(results, filename)
    # print("wrote file: %s" % filename)

    # when results are a list of dicts
    # for row in results:
    #     print(row)

    # when results are a dict
    # for key, value in results.items():
    #     print()
    #     print(key)
    #     print()
    #     print(value[0])

