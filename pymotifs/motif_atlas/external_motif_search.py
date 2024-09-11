"""
Main entry point for mapping tier 2 loops, outside the motif atlas,
to tier 1 loops, inside the motif atlas.
"""

from collections import defaultdict
from sys import path
from time import time
import os.path
import sys
import pickle
# import numpy as np
# import scipy.io as sio
# import json


# import the version of urlretrieve appropriate to the Python version
if sys.version_info[0] < 3:
    from urllib import urlretrieve as urlretrieve
    from urllib2 import urlopen
else:
    from urllib.request import urlretrieve as urlretrieve
    from urllib.request import urlopen


# this section gives us access to import things from the pythoncode folder
#wd = path[0] # current working directory "~/Dropbox/MotifAtlas"
#place = wd.rfind("\\")
#otherFolder = wd[:place] + "\\2018 FR3D Intersecting Pairs\\pythoncode"
#path.insert(1, otherFolder)

# add the search folder to the Python path
sys.path.append('search')

from fr3d_configuration import DATAPATHUNITS, DATAPATHRESULTS, DATAPATHLOOPS, SERVER
from fr3d_interactions import get_fr3d_pair_to_interaction_list
from compare_and_cluster import *  # used to be specific but copilot added things that weren't there
from test_motif_atlas_code import add_bulged_nucleotides


'''
All stacks and basepairs:
'''
bptypes = {'cWW', 'tWW', 'cWH', 'tWH', 'cWS', 'tWS', 'cHH', 'tHH', 'cHS', 'tHS',\
                   'cSS', 'tSS','cHW','tHW','cSW','tSW','cSH','tSH'}
near_bptypes = {'ntHS', 'ntSH', 'ntWH', 'ncHS', 'ncSH', 'ncWW', 'ncHH', 'ntWW', 'ntSS',\
                    'ncSS', 'ncWS', 'ntWS', 'ntHH', 'ntHW', 'ncSW', 'ncHW', 'ncWH', 'ntSW'}
stacks = {'s33','s35','s55','s53'}
near_stacks = {'ns35','ns55','ns33','ns53'}

all_stacks = stacks | near_stacks
all_bptypes = bptypes | near_bptypes
CONFLICTING_BASEPAIRS_AND_STACKS = 1
NO_CANDIDATES          = 0
SEARCH_SPACE_CONFLICT  = -1
FLANKING_BP_CONFLICT   = 9
HAIRPIN_STACK_PENALTY1 = 3
HAIRPIN_STACK_PENALTY2 = 3
BP_PENALTY             = 4
NEAR_BP_PENALTY        = 5
STACK_PENALTY          = 6
MISMATCHED_BULGE       = 10
SERVER = False    # to suppress printing list lengths
LOOP_TYPE = ""


# a temporary pretty print
def p(query):
    for key, value in query.items():
        if type(value) != defaultdict:
            print(key + ": " + str(value))
        if type(value) == defaultdict:
            print(key + ":")
            for key2, val2 in value.items():
                print("\t key {}: {}".format(key2, val2))





# just a wrapper for easy fr3d calls in testing
def wFred(Q, ifedata, listOfPairs):
    result = FR3D_search(Q, ifedata, Q['searchFiles'][0], 0)
    return(result)


def print_discrepancy_matrix(matrix, all_loops_ids):
    '''
    Helper function
    '''
    df = pd.DataFrame(matrix, index = all_loops_ids)
    with pd.option_context('display.max_rows', None,
                       'display.max_columns', None,
                       'display.precision', 3,
                       ):
        print(df)
    return()


def show_discrepancy(search_results, targets, loops_of_interest):
    dq_codes = defaultdict(list)
    discs = []
    highest_disc = 0
    highest_disc_loop = ""
    for target in targets:
        for loop in loops_of_interest:
            if loop == target:
                continue
            if (target,loop) in search_results.keys():
                dq = search_results[(target, loop)][0]['dq']
                disc = search_results[(target, loop)][0]['discrepancy']
                if dq == [-1] or dq ==[0]:
                    dq = search_results[(loop, target)][0]['dq']
                    disc = search_results[(loop, target)][0]['discrepancy']
                if len(dq) >0:
                    dq_codes[dq[0]] = loop
                else:
                    discs.append(disc)
                    if disc > highest_disc:
                        highest_disc = disc
                        highest_disc_loop = loop
        print()
        print("Testing loop:{}".format(target))
        print("Dq: {}".format(dq_codes))
        average_disc = sum(discs) / float(len(discs))
        print("Highest disc:{}, Average disc:{}".format(highest_disc, average_disc))


def startup_external_dictionaries(external_loops):
    '''
    Loop over the loop ids and extract the PDB IDs and store them in a set
    maybe also pass in "molecule type" for use in Q['requiredMoleculeType']
    '''

    all_structures = []

    for loop_id in external_loops:
            fields = loop_id.split("_")
            all_structures.append(fields[1])

    # here is where we will do file checking
    loops_and_strands = make_external_loops_and_strands(all_structures, external_loops)

    pair_to_interaction_list = defaultdict(list)
    for pdb_id in set(all_structures):
        returned_pair_to_interaction_list, returned_pair_to_crossing_number = get_fr3d_pair_to_interaction_list(pdb_id)
        returned_pair_to_interaction_list = strip_symmetry_from_pair_to_interaction(returned_pair_to_interaction_list)
        pair_to_interaction_list.update(returned_pair_to_interaction_list)



    # Loop over each loop id to identify and store the unit ids in each strand
    # loops will be a list of dictionaries, one for each loop
    # This adds a key called "strand"

    loops = strandify(loops_and_strands, all_structures) # problem
    loops = add_bulged_nucleotides(loops, pair_to_interaction_list)


    # return(loops, pair_to_interaction_list)
    return loops


def make_external_loops_and_strands(all_structures, external_loops):
    '''
    special treatment for structures from files likely external to the motif atlas
    '''

    loops_and_strands = {}

    if not os.path.exists(DATAPATHLOOPS):
        os.mkdir(DATAPATHLOOPS)

    for structure in set(all_structures):
        filename = structure + ".csv"
        # filename = structure
        pathAndFileName = os.path.join(DATAPATHLOOPS, filename)
        if not os.path.isfile(pathAndFileName):
            print("Downloading " + structure + " file...")
            urlretrieve("http://rna.bgsu.edu/rna3dhub/loops/download_with_breaks/" + structure, pathAndFileName)
            print("link: http://rna.bgsu.edu/rna3dhub/loops/download_with_breaks/" + structure)
            print("downoading TO: ", pathAndFileName)
        with open(pathAndFileName, "r") as file:
            # this next line says "for each external_loop found in this structure"
            for loop in list(filter(lambda found_loop: structure in found_loop, external_loops)):
                # Adam im expecting this to be cumbersome
                loops_and_strands[loop] = convert_csv_to_loop_data(file, loop)

    return(loops_and_strands)


def convert_csv_to_loop_data(file, loop):
    '''
    special treatment for structures not found in the motif atlas (continued)
    '''

    # go to the beginning of the file
    file.seek(0)
    internal_dict = {}

    for line in file:
        if loop in line:
            loop_nts_and_borders = line.strip("\"").replace("\"\n", "").split("\",\"")
            borders = loop_nts_and_borders[2].split(",")
            nts = loop_nts_and_borders[1].split(",")
            for index in range(len(nts)):
                internal_dict[index + 1] = (int(borders[index]), nts[index])

            return(internal_dict)


def read_results(search_results, cutoff_value = 1, require_no_dq_code = True, report_only_best_match = False):
        actual_matches = {}
        query_set = [] # starts as a list for easy appending. becomes set later

        for searched_tuple in search_results.keys():
            for match in search_results[searched_tuple]: # an element of a list, which should only ever have length one
                if(match['discrepancy'] < cutoff_value):
                    query_set.append(searched_tuple[0]) # records the query from (query, search_space)
                    if(require_no_dq_code == True):
                        if(match['dq'] == []):
                            actual_matches[searched_tuple] = search_results[searched_tuple]
                    elif(require_no_dq_code == False):
                        actual_matches[searched_tuple] = search_results[searched_tuple]

        if(report_only_best_match == False):
            for key, values in actual_matches.items():
                print()
                print(key)
                for key2, value2 in values[0].items():
                    print("\t%s: %s" % (key2, value2))

        elif(report_only_best_match == True):
            query_set = set(query_set) # removes duplicates
            best_matches = {} # will be best_match[query] = (query, search_space)
            best_rfam_matches = {}

            # begin loops to find BEST match for each query
            for query in query_set:
                best_discrepancy = 99
                best_rfam_matched_disc = 99

                for query_search_space_pair, results in actual_matches.items():
                    if(results[0]['query_id'] == query):
                        if(require_no_dq_code == True):
                            if(results[0]['dq'] == []):
                                if(results[0]['discrepancy'] < best_discrepancy):
                                    best_discrepancy = results[0]['discrepancy']
                                    best_matches[query] = query_search_space_pair
                                # if(results[0]['rfam'] not in [["disagree"], ["no rfamily"]]):
                                if('query_rfam' not in results[0]):
                                    breakpoint()

                                if(do_rfam_pools_match(results[0]['query_rfam'], results[0]['search_space_rfam'])):
                                    if('no rfam family' not in rfam_pool(results[0]['query_rfam'])):
                                        if(results[0]['discrepancy'] < best_rfam_matched_disc):
                                            best_rfam_matched_disc = results[0]['discrepancy']
                                            best_rfam_matches[query] = query_search_space_pair
                        elif(require_no_dq_code == False):
                            if(results[0]['discrepancy' < best_discrepancy]):
                                best_discrepancy = results[0]['discrepancy']
                                best_matches[query] = query_search_space_pair

            # report best matches per query
            for query, query_search_space_pair in best_matches.items():
                print()
                print(query_search_space_pair)
                for key, value in actual_matches[query_search_space_pair][0].items():
                    print("\t%s: %s" % (key, value))

            # reporting only rfam matches for debugging
            for query, query_search_space_pair in best_rfam_matches.items():
                print()
                print(query_search_space_pair)
                for key, value in actual_matches[query_search_space_pair][0].items():
                    print("\t%s: %s" % (key, value))


def load_previous_search_results(query_pdb_id, path = DATAPATHRESULTS, loop_type = "IL"):

    results = {}
    if not os.path.exists(path):
        os.mkdir(path)
        return(results)
    #file_name_and_path = path + "/search_results/" + structure_name + ".pickle"
    file_name = loop_type + "_" + query_pdb_id + "_search_results.pickle"
    file_path = os.path.join(path,file_name)

    if os.path.exists(file_path):
        if sys.version_info[0] < 3:
            results = pickle.load(open(file_path, "rb"))
        else:
            results = pickle.load(open(file_path, "rb"), encoding = 'latin1')

        for (query,search_space),value in results.items(): #Load No-match result as [{'dq': 0, 'discrepancy': 99}] instead of an integer
            if isinstance(value,int):
                results[(query,search_space)] = [{'dq': [value], 'discrepancy': 99}]

        return(results)
    return(results)


def write_results_manual(search_results, cutoff_value = 1, require_no_dq_code = True, filename = "adam_dont.csv"):
    import csv

    actual_matches = {}
    query_set = [] # starts as a list for easy appending. becomes set later

    # removes results below a cutoff value, and if there are dq codes
    for searched_tuple in search_results.keys():
        for match in search_results[searched_tuple]: # an element of a list, which should only ever have length one
            if(match['discrepancy'] < cutoff_value):
                query_set.append(searched_tuple[0]) # records the query from (query, search_space)
                if(require_no_dq_code == True):
                    if(match['dq'] == []):
                        actual_matches[searched_tuple] = search_results[searched_tuple]
                elif(require_no_dq_code == False):
                    actual_matches[searched_tuple] = search_results[searched_tuple]


    query_set = set(query_set) # removes duplicates
    best_matches = {} # will be best_match[query] = (query, search_space)
    best_rfam_matches = {}


    # begin loops to find BEST match for each query
    for query in query_set:
        best_discrepancy = 99
        best_rfam_matched_disc = 99

        for query_search_space_pair, results in actual_matches.items():
            if(results[0]['query_id'] == query):
                if(require_no_dq_code == True):
                    if(results[0]['dq'] == []):

                        if(results[0]['discrepancy'] < best_discrepancy):
                            best_discrepancy = results[0]['discrepancy']
                            best_matches[query] = query_search_space_pair

                        if(do_rfam_pools_match(results[0]['query_rfam'], results[0]['search_space_rfam'])):
                            if('no rfam family' not in rfam_pool(results[0]['query_rfam'])):
                                if(results[0]['discrepancy'] < best_rfam_matched_disc):
                                    best_rfam_matched_disc = results[0]['discrepancy']
                                    best_rfam_matches[query] = query_search_space_pair
                elif(require_no_dq_code == False):
                    if(results[0]['discrepancy' < best_discrepancy]):
                        best_discrepancy = results[0]['discrepancy']
                        best_matches[query] = query_search_space_pair

    with open(filename, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=',')
        # old header = ["Query", "Search Space", "Rfams Match", "Discrepancy", "q_rfam", "ss_rfam"]
        # new header indices 0 = Q, 1 = SS, 2 = Disc, 3 = Rf match, 4 = EC match, 5 = Q RF, 6 = SS RF, 7 = Q EC, 8 = SS EC
        header = ["Query", "Search Space", "Discrepancy", "Rfams Match", "Equiv Class Match", "q_rfam", "ss_rfam", "q_equiv", "ss_equiv"]
        csvwriter.writerow(header)

        # report best geometric matches per query
        for query, query_search_space_pair in best_matches.items():
            # breakpoint()

            # the "missing" processing in the next 2 blocks should be moved to INSIDE the functions they call
            # if(do_rfam_pools_match(results[0]['query_rfam'], results[0]['search_space_rfam'])):
            if(do_rfam_pools_match(search_results[query_search_space_pair][0]['query_rfam'], search_results[query_search_space_pair][0]['search_space_rfam'])):
                # if('no rfam family' not in rfam_pool(results[0]['query_rfam'])):
                if('no rfam family' not in rfam_pool(search_results[query_search_space_pair][0]['query_rfam'])):
                    homologous = True
                else:
                    homologous = False
            else:
                homologous = False

            if(do_equiv_pools_match(search_results[query_search_space_pair][0]['query_rfam'], search_results[query_search_space_pair][0]['search_space_rfam'])):
                if("no equivalence class" not in equiv_pool(search_results[query_search_space_pair][0]['query_rfam'])):
                    equivalence_match = True
                else:
                    equivalence_match = False
            else:
                equivalence_match = False

            search_space = query_search_space_pair[1]
            discrepancy = actual_matches[query_search_space_pair][0]["discrepancy"]
            query_rfam = organize_family_per_chain(actual_matches[query_search_space_pair][0]["query_rfam"])
            search_space_rfam = organize_family_per_chain(actual_matches[query_search_space_pair][0]["search_space_rfam"])
            query_equiv_class = organize_equiv_per_chain(actual_matches[query_search_space_pair][0]["query_rfam"])
            search_space_equiv_class = organize_equiv_per_chain(actual_matches[query_search_space_pair][0]["search_space_rfam"])


            # old row = [query, search_space, discrepancy, homologous, actual_matches[query_search_space_pair][0]["query_rfam"], actual_matches[query_search_space_pair][0]["search_space_rfam"]]
            row = [query, search_space, discrepancy, homologous, equivalence_match, query_rfam, search_space_rfam, query_equiv_class, search_space_equiv_class]
            csvwriter.writerow(row)

        csvwriter.writerow([""]) # empty row to see difference

        # reporting only rfam matches for debugging
        for query, query_search_space_pair in best_rfam_matches.items():
            homologous = True # this stems from the requirement to get into this dict

            if(do_equiv_pools_match(results[query_search_space_pair]['query_rfam'], results[query_search_space_pair]['search_space_rfam'])):
                if("no equivalence class" not in equiv_pool(results[query_search_space_pair]['query_rfam'])):
                    equivalence_match = True
                else:
                    equivalence_match = False
            else:
                equivalence_match = False

            search_space = query_search_space_pair[1]
            discrepancy = actual_matches[query_search_space_pair][0]["discrepancy"]
            query_rfam = organize_family_per_chain(actual_matches[query_search_space_pair][0]["query_rfam"])
            search_space_rfam = organize_family_per_chain(actual_matches[query_search_space_pair][0]["search_space_rfam"])
            query_equiv_class = organize_equiv_per_chain(actual_matches[query_search_space_pair][0]["query_rfam"])
            search_space_equiv_class = organize_equiv_per_chain(actual_matches[query_search_space_pair][0]["search_space_rfam"])

            # old row = [query, query_search_space_pair[1], homologous, actual_matches[query_search_space_pair][0]["discrepancy"], actual_matches[query_search_space_pair][0]["query_rfam"], actual_matches[query_search_space_pair][0]["search_space_rfam"]]
            row = [query, search_space, discrepancy, homologous, equivalence_match, query_rfam, search_space_rfam, query_equiv_class, search_space_equiv_class]
            csvwriter.writerow(row)


def filter_to_only_best_chains(current_chains = []):
    '''
    narrowing down our chains of interest by composite quality score.
    If 2 of the chains here are in the same equivalence class, the one with
    less quality should be omitted
    '''

    chain_to_equiv_class = read_equiv_class_csv_into_dict()
    equiv_class_to_chains = {}
    best_chains = []
    best_candidate_chain = ''

    for chain in current_chains:
        if chain_to_equiv_class[chain]['equivalence class'] not in equiv_class_to_chains:
            equiv_class_to_chains[chain_to_equiv_class[chain]['equivalence class']] = [chain]
        else:
            equiv_class_to_chains[chain_to_equiv_class[chain]['equivalence class']].append(chain)

    for equiv_class, list_of_chains in equiv_class_to_chains.items():
        if len(list_of_chains) == 1:
            best_chains.append(list_of_chains[0])
        else:
            lowest_rank = 50000000 # lower is better, this sentinel value is fragile and needed raise twice
            for chain in list_of_chains:
                # find lowest rank to make best
                if chain_to_equiv_class[chain]['quality rank'] < lowest_rank:
                    lowest_rank = chain_to_equiv_class[chain]['quality rank']
                    best_candidate_chain = chain
            best_chains.append(best_candidate_chain)

    return(best_chains)


def make_list_of_all_loops_in(structure, number_of_strands = 0, max_size = None):
    '''
    fetches/opens a file of all the loops in a given structure, then does a filtered read
    based on the number of strands that we care about, to return only the loops with that
    many strands participating (IE: 2 strands = internal loops, 1 strand = hairpin loops)

    structure is also allowed to be a chain id or IFE
    example: structure = 6XU8|1|A5+6XU8|1|A8
    works =]
    '''

    loops = []
    structures = structure.split('+')

    chain_names = set([chain.split("|")[2] for chain in structures])

    if not os.path.exists(DATAPATHUNITS):
        os.mkdir(DATAPATHUNITS) # if there's an error here, find the value of DATAPATHUNITS and manually make it

    if number_of_strands == 1:
        prefix = "HL_"
    elif number_of_strands == 2:
        prefix = "IL_"
    elif number_of_strands > 2:
        prefix = "J%d_" % number_of_strands
    elif number_of_strands == 0:
        prefix = ""

    pdb_id = structure.split("|")[0]
    filename = pdb_id + ".csv"

    pathAndFileName = os.path.join(DATAPATHLOOPS, filename)
    if not os.path.isfile(pathAndFileName):
        print("Downloading " + pdb_id + " file...")
        urlretrieve("http://rna.bgsu.edu/rna3dhub/loops/download_with_breaks/" + pdb_id, pathAndFileName)

    with open(pathAndFileName, "r") as file:
        # add a clause here where if no loops found
        # tries to redownload
        for line in file:
            fields = line.split(",")
            if len(fields) < 2:
                continue
            loop_id = fields[0].strip('"')
            if number_of_strands == 0 or prefix in loop_id:
                unit_id = fields[1].strip('"')
                c_fields = unit_id.split("|")
                if len(c_fields) >= 3:
                    chain_name = c_fields[2]
                    if chain_name in chain_names:
                        if not max_size or len(fields) - 1 <= 2*max_size:
                            loops.append(loop_id)

    return list(set(loops))


# IL_6ZMI_0048 has missing nts, IL_6ZMI_195 has incomplete nts. http://rna.bgsu.edu/rna3dhub/pdb/6ZMI/motifs
def main(loop_ids_of_interest = ["IL_6ZMI_008", "IL_6ZMI_007", "IL_6ZMI_048", "IL_6ZMI_195", "IL_4YAZ_008"], save_path = os.path.join(DATAPATHRESULTS,"for_annotation"), load_previous_result = True):
    '''
    it will try to align all of the loop_ids_of_interest to the
    loops in the motif atlas
    '''

    from search import FR3D_search, matrix_discrepancy_cutoff
    from file_reading import readNAPairsFile
    from file_reading import readProteinPositionsFile
    from myTimer import myTimer
    # using these to get other args to feed FR3D_search()
    # from pair_processing import get_pairlist # this got moved into search.py
    from query_processing import calculateQueryConstraints
    from query_processing import retrieveQueryInformation
    from chain_to_rfam_family import do_rfam_pools_match, do_equiv_pools_match, rfam_pool, equiv_pool, read_equiv_class_csv_into_dict, organize_equiv_per_chain, organize_family_per_chain
    # from query_processing import emptyInteractionMatrix
    # a program on the server can generate the next line
    # That lists all the loops in a motif atlas and the unit ids in each strand.  I think.
    # Keys are loop ids, then positions and maybe an indication of what strand they are in
    # (1, '6DVK|1|H|G|28') means that G28 is on the border of a single-stranded region
    # Then A29,
    # http://rna.bgsu.edu/rna3dhub/loops/view/IL_6DVK_003

    from IL_3_57_loops_and_strands import loops_and_strands

    timer_data = myTimer('Setup')

    # tier 1 means "from inside the motif atlas"
    # tier 2 means "from outside of the motif atlas"

    # these are setup up in a parallel manner, for the scenario in which we reverse the search
    tier_1_queries = {}
    tier_1_search_spaces = {}
    tier_1_flanking_bp_queries = {}

    tier_2_queries = {}
    tier_2_search_spaces = {}
    tier_2_flanking_bp_queries = {}

    # THIS IS VERY TEMPORARY #################################################################
    # PROOF OF CONCEPT FOR STRUCTURE VS STRUCTURE:
    # new_filtered = {}
    # for loop in loops_and_strands:
    #     if("7RQB" in loop):
    #         new_filtered[loop] = loops_and_strands[loop]

    #loops_and_strands = new_filtered
    ###########################################################################################

    # timer_data = myTimer('dicts')
    tier_1_loops, trash = startup_list_of_dictionaries(loops_and_strands)
    # tier_1_loops, trash = startup_list_of_dictionaries(new_filtered)
    tier_2_loops = startup_external_dictionaries(loop_ids_of_interest)
    # breakpoint()

    # OLD, delete once flanking queries are made inside "create_q_and_ss()"
    # # timer_data = myTimer('queries and search spaces')
    # tier_1_queries, tier_1_search_spaces = create_queries_and_search_spaces(tier_1_loops)
    # tier_2_queries, tier_2_search_spaces = create_queries_and_search_spaces(tier_2_loops, loop_ids_of_interest)
    # # timer_data = myTimer('flanking bp queries?')
    # tier_1_flanking_bp_queries = create_flanking_bp_queries(tier_1_loops)
    # tier_2_flanking_bp_queries = create_flanking_bp_queries(tier_2_loops, loop_ids_of_interest)

    # timer_data = myTimer('queries and search spaces')
    tier_1_queries, tier_1_search_spaces, tier_1_flanking_bp_queries = create_queries_and_search_spaces(tier_1_loops)
    tier_2_queries, tier_2_search_spaces, tier_2_flanking_bp_queries = create_queries_and_search_spaces(tier_2_loops, loop_ids_of_interest)
    breakpoint()
    set_loop_type(tier_2_queries) # this sets a global variable. Bad practice

    # this is "comparing and aligning" the tier 2 loops as queries, in the tier 1 search spaces
    search_results = all_against_all_searches(tier_2_queries, tier_1_search_spaces, tier_2_flanking_bp_queries, load_previous_result = load_previous_result, reversed_search = False, save_path = save_path)

    # print("Finished searches, beginning output")
    timer_data = myTimer("Output")

    print(myTimer('summary'))

    read_results(search_results, report_only_best_match = True)

    return(search_results)





# if __name__ == '__main__':
    # human_ribosome = make_list_of_all_loops_in(structure = "6ZMI", number_of_strands = 2)
    # results = main(loop_ids_of_interest = human_ribosome)
    # main()