"""
the main purpose of this script is to find annotations for
loops, given a chain
by Adam Coger, acoger@bgsu.edu 04/24/2023
and Craig Zirbel, zirbel@bgsu.edu 01/25/2025
"""

import copy
import json
import os
import pickle
import sys
import time
from urllib.request import urlopen

from fr3d.search.search import FR3D_search, lookUpInteractions

# from pymotifs.motif_atlas.clustering_utilities import load_motif_atlas_release
from pymotifs.motif_atlas.compare_and_cluster import create_queries_and_search_spaces
from pymotifs.motif_atlas.compare_and_cluster import list_all_chains_in
from pymotifs.motif_atlas.compare_and_cluster import nt_to_loop_mapping
from pymotifs.motif_atlas.compare_and_cluster import filter_on_unmatched_nucleotides
from pymotifs.motif_atlas.compare_and_cluster import filter_out_conflicting_basepairs_and_stacks
from pymotifs.motif_atlas.compare_and_cluster import keep_lowest_discrepancy_candidate
from pymotifs.motif_atlas.compare_and_cluster import all_against_all_searches
from pymotifs.motif_atlas.compare_and_cluster import load_previous_search_results_one_pdb
from pymotifs.motif_atlas.external_motif_search import make_list_of_all_loops_in
from pymotifs.motif_atlas.external_motif_search import startup_external_dictionaries
from pymotifs.motif_atlas.myTimer import myTimer

# from pymotifs.motif_atlas.compare_motifs_of_two_structures import write_output_csv

from pymotifs.motif_atlas.discrepancy_flip import matrix_discrepancy_cutoff_flip as matrix_discrepancy_cutoff
# from fr3d.geometry.discrepancy import matrix_discrepancy_cutoff

DATAPATHATLAS = '/usr/local/pipeline/hub-core/MotifAtlas/tier2/'
DATAPATHRESULTS = '/usr/local/pipeline/hub-core/MotifAtlas/tier2/search_results'
DATAPATHALIGNMENTS = '/usr/local/pipeline/hub-core/MotifAtlas/tier2/alignments'
DATAPATHMOTIFDATA = os.path.join(DATAPATHATLAS,'motif_data')

SERVER = True

homologous_discrepancy_cutoff = 0.6
geometric_discrepancy_cutoff = 0.4


def save_motif_atlas_already_set_up(release, atlas_queries, atlas_search_spaces, atlas_flanking_bp_queries):
    """
    this saves some complicated dictionaries to avoid constant reprocessing
    """
    if not os.path.exists(DATAPATHMOTIFDATA) and not SERVER:
        os.mkdir(DATAPATHMOTIFDATA)

    pickle.dump(obj = atlas_queries, file = open(os.path.join(DATAPATHMOTIFDATA,
        "%s_atlas_queries.pickle" % release), "wb"), protocol = 4)
    pickle.dump(obj = atlas_search_spaces, file = open(os.path.join(DATAPATHMOTIFDATA,
        "%s_atlas_search_spaces.pickle" % release), "wb"), protocol = 4)
    pickle.dump(obj = atlas_flanking_bp_queries, file = open(os.path.join(DATAPATHMOTIFDATA,
        "%s_atlas_flanking_bp_queries.pickle" % release), "wb"), protocol = 4)
    return


def load_motif_atlas_already_set_up(release):
    """
    this is the counterpart to the above function "save_motif_atlas_already_set_up"
    """
    atlas_queries = pickle.load(open(os.path.join(DATAPATHMOTIFDATA,
        "%s_atlas_queries.pickle" % release), "rb"))
    atlas_search_spaces = pickle.load(open(os.path.join(DATAPATHMOTIFDATA,
        "%s_atlas_search_spaces.pickle" % release), "rb"))
    atlas_flanking_bp_queries = pickle.load(open(os.path.join(DATAPATHMOTIFDATA,
        "%s_atlas_flanking_bp_queries.pickle" % release), "rb"))

    return atlas_queries, atlas_search_spaces, atlas_flanking_bp_queries


# def get_current_motif_atlas_ids(types = ["HL", "IL", "J3"], release = "current"):
#     """
#     Download json object and parse to get loop ids
#     """

#     loop_ids = []
#     motif_list = []

#     for loop_type in types:
#         print("Loading motif atlas release %s for %s" % (release,loop_type))
#         motif_list.extend(load_motif_atlas_release(loop_type = loop_type, release = release, path = DATAPATHATLAS))

#     for motif in motif_list:
#         loop_ids.extend(motif['alignment'].keys())

#     return loop_ids


def get_motif_atlas_release_loop_ids(loop_type = 'HL', release = 'current'):
    """
    Download motif atlas release
    Return release id and loop ids
    """

    url = "https://rna.bgsu.edu/rna3dhub/motifs/release/%s/%s/json" % (loop_type,release)
    response = urlopen(url)

    print('Downloading %s' % url)

    # Get the filename from the Content-Disposition header
    content_disposition = response.headers.get('Content-Disposition')
    if content_disposition and 'filename=' in content_disposition:
        filename = content_disposition.split('filename=')[-1].strip(' "')
        release = filename.replace(".json","").split("_")[1]

    # Read the content
    content = response.read().decode('utf-8')

    try:
        motif_groups = json.loads(content)
    except json.JSONDecodeError:
        motif_groups = []

    loop_ids = []
    for motif_group in motif_groups:
        if 'alignment' in motif_group:
            if type(motif_group['alignment']) == dict:
                loop_ids.extend(motif_group['alignment'].keys())

    return release, loop_ids


def update_motif_atlas_files(release='current'):
    """
    Check that we have created the specified motif atlas release files to run against
    """

    release, motif_atlas_loop_ids = get_motif_atlas_release_loop_ids('HL',release)

    # check if this release has already been downloaded and processed
    files_already_exist = bool(os.path.isfile(os.path.join(DATAPATHMOTIFDATA, "%s_atlas_queries.pickle" % release)) and
        os.path.isfile(os.path.join(DATAPATHMOTIFDATA, "%s_atlas_search_spaces.pickle" % release)) and
        os.path.isfile(os.path.join(DATAPATHMOTIFDATA, "%s_atlas_flanking_bp_queries.pickle" % release)))

    if not files_already_exist:
        print("Creating motif atlas files for release %s" % release)
        additional_loop_types = ["IL", "J3"]    # recently added J3 but might cause an error
        for loop_type in additional_loop_types:
            release, loop_ids = get_motif_atlas_release_loop_ids(loop_type,release)
            motif_atlas_loop_ids.extend(loop_ids)

        print('Getting loop data')
        loops_of_atlas = startup_external_dictionaries(motif_atlas_loop_ids)
        print('Making queries and search spaces')
        atlas_queries, atlas_search_spaces, atlas_flanking_bp_queries, loops_of_atlas = create_queries_and_search_spaces(loops_of_atlas,geometric_discrepancy_cutoff)
        print('Saving queries and search spaces')
        save_motif_atlas_already_set_up(release, atlas_queries, atlas_search_spaces, atlas_flanking_bp_queries)

    return release


def calculate_discrepancy_from_alignment(Q, ifedata, timerData):

    # compute the discrepancy between the nucleotides in Q and the nucleotides in ifedata
    # the two sets of nucleotides are assumed to be in the same order, aligned

    # It seems like the plan was to pass in a permutation to tell how to align them,
    # but that is not happening now.

    # perm comes from listOfPairs from ifedata
    # perm = range(0, 48 + 1) # positions of non-bulged and aligned nucleotides
    perm = range(0, len(ifedata['centers'])) # positions of non-bulged and aligned nucleotides
    # perm is saying which nucleotides to use
    # might want to pass in perm. which is a list of indices to get compared
    # that list perm would say which indices of nts from ifedata will be used

    positions = list(range(0, Q['numPositions'])) # set for straight alignment

    querycenters = [Q["centers"][i] for i in positions]
    queryrotations = [Q["rotations"][i] for i in positions]
    # querycenters = [Q["centers"][i] for i in perm]
    # queryrotations = [Q["rotations"][i] for i in perm]

    units = ifedata['units']
    index_to_id = ifedata['index_to_id']
    pairToInteractions = ifedata['pairToInteractions']
    pairToCrossingNumber = ifedata['pairToCrossingNumber']
    timerData = myTimer("Discrepancy from query")

    ifedata_centers = []
    for i in range(0, Q['numPositions']):
        ifedata_centers.append(units[i]["centers"])
    ifedata_rotations = []
    for i in range(0, Q['numPositions']):
        ifedata_rotations.append(units[i]["rotations"])

    d = matrix_discrepancy_cutoff(querycenters, queryrotations, ifedata_centers,
        ifedata_rotations, Q["discrepancy"])

    candidates = []

    if d is not None and d <= Q["discrepancy"]:
        newcandidate = {}
        indices = positions
        newcandidate['indices'] = indices
        newcandidate['unitids'] = [index_to_id[index] for index in indices]
        newcandidate['chainindices'] = [units[index]["chainindex"] for index in indices]
        newcandidate['centers'] = [units[index]["centers"] for index in indices]
        newcandidate['rotations'] = [units[index]["rotations"] for index in indices]
        newcandidate['discrepancy'] = d
        newcandidate['interactions'] = lookUpInteractions(Q,indices,
            pairToInteractions, pairToCrossingNumber, units)
        candidates.append(newcandidate)

    return Q, candidates, timerData


def run_aligned_searches_first(queries = {}, search_spaces = {}, structure_a_best_chains = False, flanking_bp_queries = {}, return_all_queries = False, force_alignment = False):
    """
    Use a chain to chain alignment to suggest homologous motif matches,
    check the geometric discrepancy.
    Matches are likely because of homology.
    """
    # queries = structure B, search_spaces = structure A

    # Aligned loops can have discrepancy up to this number
    homologous_discrepancy_cutoff = 0.6

    aligned_search_pairs = {}
    # queries_to_still_run = {}
    search_spaces_without_matches = {}
    search_results = {}
    failed_aligned_search_pairs = {}
    Q = {}

    # could be a way to focus down list of chains in a
    query_chains = list_all_chains_in(queries)

    # best chains stuff should be done before this function is called
    if structure_a_best_chains:
        search_space_chains = structure_a_best_chains
    else:
        search_space_chains = list_all_chains_in(search_spaces)

    nt_to_loop_of_ss = nt_to_loop_mapping(search_spaces)
    nt_to_loop_of_q = nt_to_loop_mapping(queries)

    # filter chain list of search space to "best chains"
    # some_new_function_to_do_that()

    # if not os.path.exists(DATAPATHALIGNMENTS):
    #     os.mkdir(DATAPATHALIGNMENTS)
    # files_of_dir = os.listdir(DATAPATHALIGNMENTS)
    # txt_files = [file for file in files_of_dir if ".txt" in file]

    for ss_chain in search_space_chains:
        for q_chain in query_chains:
            # map the loops of these chains will need to return some alignment data
            loop_to_aligned_loop = map_the_loops_of_these_chains(queries, search_spaces, q_chain, ss_chain, nt_to_loop_of_ss)
            # print("annotations_for_chain_loops: loop_to_aligned_loop:")
            # print(loop_to_aligned_loop)
            aligned_search_pairs.update(loop_to_aligned_loop)
    timer_data = myTimer("Running Aligned Searches")

    for query_id, search_space_id in sorted(aligned_search_pairs.items()):
        print("Doing aligned search of %s in %s" % (query_id, search_space_id))
        #print("Query:",queries[query_id]['Q']['unitID'])

        # this is later than this should be set; set it when the query is made, to have
        # tighter constraints on the whole search
        # queries[query_id]['Q']['discrepancy'] = homologous_discrepancy_cutoff

        search_space_length = len(search_spaces[search_space_id]['ifedata']['centers'])
        candidates = []

        query_unit_id = queries[query_id]['Q']['unit_id']

        # if we can align directly from query (or query with bulges) to search space, try that
        if queries[query_id]['Q']['numPositions'] == search_space_length:
            # non-bulged nucleotides from the query aligned to all search space nucleotides
            print('Using non-bulged query nucleotides')
            Q, candidates, elapsed_time = calculate_discrepancy_from_alignment(Q = queries[query_id]['Q'],
                ifedata = search_spaces[search_space_id]['ifedata'],
                timerData = timer_data)
        elif len(queries[query_id]['Q']['fullUnits']) == search_space_length:
            # all nucleotides from the query, including bulges, aligned to all search space nucleotides

            # make a copy of the query so we can use all nucleotides including bulges
            mockQ = copy.deepcopy(queries[query_id]['Q'])
            mockQ['numPositions'] = len(mockQ['fullUnits'])
            mockQ['centers'] = queries[query_id]['ifedata']['centers']
            mockQ['rotations'] = queries[query_id]['ifedata']['rotations']

            Q, candidates, elapsed_time = calculate_discrepancy_from_alignment(Q = mockQ,
                ifedata = search_spaces[search_space_id]['ifedata'],
                timerData = timer_data)

            query_unit_id = queries[query_id]['Q']['fullUnits']
            print('Using all query nucleotides for %s with %s' % (query_id,search_space_id))

        if not candidates:
            # if no alignment, or alignment did not give a low discrepancy, do a regular fr3d search
            if search_space_length < 30:
                Q, candidates, elapsed_time = FR3D_search(Q = queries[query_id]['Q'],
                                                ifedata = search_spaces[search_space_id]['ifedata'],
                                                ifename = search_space_id,
                                                timerData = timer_data)

        if candidates:
            temp_result = filter_out_conflicting_basepairs_and_stacks(candidates, queries[query_id], search_spaces[search_space_id], query_unit_id)
            if temp_result:
                temp_result = filter_on_unmatched_nucleotides(temp_result,search_spaces[search_space_id],queries[query_id])
                temp_result = keep_lowest_discrepancy_candidate(temp_result)

                if temp_result:
                    if temp_result[0]["discrepancy"] <= homologous_discrepancy_cutoff:
                        temp_result[0]['match_type'] = "homologous"
                        search_results[(query_id,search_space_id)] = temp_result
                        print('Aligned loops %s and %s match' % (query_id,search_space_id))
                        continue

        failed_aligned_search_pairs[query_id] = aligned_search_pairs[query_id]

    # ss dict minus aligned_search_pairs keys
    for loop_id, ss_info in search_spaces.items():
        if loop_id not in aligned_search_pairs.values():
            search_spaces_without_matches[loop_id] = ss_info
        if loop_id in failed_aligned_search_pairs.values():
            search_spaces_without_matches[loop_id] = ss_info

    return search_results, search_spaces_without_matches


def get_nt_alignment_between_chains(chain_of_a, chain_of_b):
    """
    uses one of our APIs like: http://rna.bgsu.edu/correspondence/pairwise_structure?chain1=4V88|1|A6&chain2=5TBW|1|sR
    which responds with the nucleotide to nucleotide alignment of two chains.
    From there that mapping is used to align motifs from A to motifs of B, and return that mapping
    """

    ntB_to_ntA_alignments = {}

    possible_file_name = "%s_to_%s_alignment.pickle" % (chain_of_a.replace("|", "-"), chain_of_b.replace("|", "-"))
    # path_and_file = os.path.join(".\\alignments", possible_file_name) # commented out 06/18/2023
    path_and_file = os.path.join(DATAPATHALIGNMENTS, possible_file_name)

    # if folder needs writing
    if not os.path.exists(DATAPATHALIGNMENTS) and not SERVER:
        os.mkdir(DATAPATHALIGNMENTS)

    # if we have the alignment. load it and move on
    if os.path.exists(path_and_file):
        ntB_to_ntA_alignments = pickle.load(open(path_and_file, 'rb'))
        return ntB_to_ntA_alignments

    url = "https://rna.bgsu.edu/correspondence/align_chains?chains=" + chain_of_b + "," + chain_of_a

    try:
        reply = urlopen(url).read()
        print("success reading %s" % url)
        reply = reply.decode("ascii")
    except:
        reply = ""
        print("failure reading %s" % url)

    for index, line in enumerate(reply.split("\n")):
        fields = line.split()
        if len(fields) == 2:
            ntB, ntA = fields[0], fields[1]
            if ntB == "NULL" or ntA == "NULL":
                pass
            else:
                # there ARE possible overwrites here =/
                ntB_to_ntA_alignments[ntB] = ntA

    # save this as a pickle file for next time
    with open(path_and_file, 'wb') as f:
        pickle.dump(ntB_to_ntA_alignments, f)

    return ntB_to_ntA_alignments


def map_the_loops_of_these_chains(queries, search_spaces, q_chains, ss_chains, nt_to_loop_of_ss):
    """
    this is mapping loops to loops from two different chains in the same Rfam family
    `ss` is short for search space, `q` is short for query
    """

    loop_to_aligned_loop = {}
    ntB_to_ntA_alignment = {}

    # ntB_to_ntA_alignments[structure a nt unit id] = structure b nt unit id
    # good to know, but that seems backward

    for q_chain in q_chains.split("+"):
        for ss_chain in ss_chains.split("+"):
            ntB_to_ntA_alignment.update(get_nt_alignment_between_chains(chain_of_a = ss_chain, chain_of_b = q_chain))

    if len(ntB_to_ntA_alignment) == 0:
        return loop_to_aligned_loop

    for loop_name_in_q, q_loop_data in queries.items():
        query_loop_type = loop_name_in_q.split("_")[0]
        # query_loop_data['loop_info']['strand']] <- this is a list of lists, ex: [["6DVK|1|H|G|28", "6DVK|1|H|A|29"], ["6DVK|1|H|G|61", "6DVK|1|H|C|62"]]
        # when a nt from search_space is aligned to a nt in a query loop, that loop's name is appended

        # this next line is a list comprehension, to flatten a list of lists into a regular 1D list
        nts_of_q_loop = [nt for strand in q_loop_data['loop_info']['strand'] for nt in strand]
        ss_loops_with_aligned_nts = []

        for nt in nts_of_q_loop:
            # new
            ss_loop = nt_to_loop_of_ss.get(ntB_to_ntA_alignment.get(nt))
            if ss_loop:
                if ss_loop.split("_")[0] == query_loop_type:
                    ss_loops_with_aligned_nts.append(ss_loop)
            # old
            # if nt in ntB_to_ntA_alignment:
            #   if ntB_to_ntA_alignment[nt] in nt_to_loop_of_ss:
            #       ss_loops_with_aligned_nts.append(nt_to_loop_of_ss[ntB_to_ntA_alignment[nt]])

        # next line find the most common match
        # max(set(query_loop_matches), key = query_loop_matches.count)
        if ss_loops_with_aligned_nts:
            most_common_ss_loop = max(set(ss_loops_with_aligned_nts), key = ss_loops_with_aligned_nts.count)

            # is the most common match over HALF of the nts?
            if ss_loops_with_aligned_nts.count(most_common_ss_loop)/len(ss_loops_with_aligned_nts) > .5:
                # add it to the loop to loop mapping for return
                loop_to_aligned_loop[loop_name_in_q] = most_common_ss_loop

    return loop_to_aligned_loop


def main(chain_id, chains_for_comparing_against = [], finished_loops = [], max_size = None, release = None):
    # chain_id is allowed to be an IFE

    # if someone input a string as the second variable, convert to a list
    if isinstance(chains_for_comparing_against, str):
        chains_for_comparing_against = [chains_for_comparing_against]

    # get loop ids from chain_id which is to be searched; preparing chain loops for FR3D / AvA searches
    loop_ids_of_chain = make_list_of_all_loops_in(structure = chain_id, max_size = max_size)

    print('Found %d loops in chain %s' % (len(loop_ids_of_chain), chain_id))

    # set aside the loops that are already mapped, if any
    loop_ids_of_chain = sorted(set(loop_ids_of_chain) - set(finished_loops))

    print('Still need to map %d loops in chain %s' % (len(loop_ids_of_chain), chain_id))

    if len(loop_ids_of_chain) == 0:
        return []

    # print("Loops in chain %s are %s" % (chain_id,sorted(loop_ids_of_chain)))

    # look up positions, border, identify strands, pass back full information
    loops_of_chain = startup_external_dictionaries(loop_ids_of_chain)

    # make the FR3D queries and search spaces for loops in the new chain to be annotated
    # however, we only use the search spaces here
    chain_as_queries, chain_as_search_spaces, chain_flanking_bp_queries, loops_of_chain = create_queries_and_search_spaces(loops_of_chain)

    if len(chain_as_search_spaces) == 0:
        return []

    problems = ["HL_6ZMI_059", "J3_6ZMI_003"]
    for problem in problems:
        if problem in chain_as_queries:
            del chain_as_queries[problem]
            del chain_as_search_spaces[problem]
            del chain_flanking_bp_queries[problem]

    print("Loading pre-computed motif atlas files for release %s" % release)
    atlas_queries, atlas_search_spaces, atlas_flanking_bp_queries = load_motif_atlas_already_set_up(release)

    # loop over chains to compare against
    double_filtered_results = {}
    aligned_results = {}
    unmatched_search_spaces = {} # this is the overlap of leftover_search_spaces from ALL annotated chains for comparing against
    search_results = {} # populated from motif atlas geometric searches

    for annotated_chain_for_comparing in chains_for_comparing_against:
        # print("in comparison loop. checking %s against %s" % (chain_id, annotated_chain_for_comparing))
        # gather loops from this annotated chain
        loop_ids_of_annotated_chain = make_list_of_all_loops_in(structure = annotated_chain_for_comparing)
        loops_of_annotated_chain = startup_external_dictionaries(loop_ids_of_annotated_chain)
        annotated_chain_as_queries, annotated_chain_as_search_spaces, annotated_chain_flanking_bp_queries, loops_of_annotated_chain = create_queries_and_search_spaces(loops_of_annotated_chain,homologous_discrepancy_cutoff)

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
                save_path = DATAPATHRESULTS)
            search_results.update(loop_type_results)

        if len(loop_type_search_spaces) < 2:
            # if there are no search spaces fed to all_against_all, it will not load prev_results
            print("in force load for %s loops" % prefix)
            results = load_previous_search_results_one_pdb(loop_type = prefix, pdb_id = chain_id.split("|")[0], path = DATAPATHRESULTS)
            for pair, info in results:
                if pair[0] in atlas_queries and pair[1] in chain_as_search_spaces:
                    if not isinstance(results[(query_id, search_space_id)],int): # says "if not disqualified"
                        # print("added a match!")
                        search_results[(query_id,search_space_id)] = results[(query_id,search_space_id)]
        # print(" %s : %d" % (prefix, len(loop_type_search_spaces)))

    # cutting out disqualified and high discrepancy results
    filtered_results = {}
    geometric_discrepancy_cutoff = 0.4
    for pair, info in search_results.items():
        if pair[1] in chain_as_search_spaces:
            if info[0]['dq'] == []:
                if info[0]['discrepancy'] <= geometric_discrepancy_cutoff:
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

    return results_in_list_of_dicts_format



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

