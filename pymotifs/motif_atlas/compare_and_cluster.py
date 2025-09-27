"""
The core of the motif atlas calculations.
Given a list of HL, IL, J3, or other loops,
do all-against-all FR3D searches to find all possible matches.
Then, cluster the loops based on the matches.
Return a set of motif groups.
"""

from collections import defaultdict, OrderedDict
import gzip
import networkx as nx
import numpy as np
import os.path
import pickle
import shutil
import sys
from sys import maxsize
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from time import time

# import the version of urlretrieve appropriate to the Python version
if sys.version_info[0] < 3:
    from urllib import urlretrieve as urlretrieve
else:
    from urllib.request import urlretrieve as urlretrieve

from pymotifs.motif_atlas.fr3d_interactions import get_fr3d_pair_to_interaction_list
from pymotifs.motif_atlas.myTimer import myTimer
from pymotifs.motif_atlas.clustering_utilities import writeCSVOutput
from pymotifs.motif_atlas.clustering_utilities import get_matrix_for_consensus_interactions
from pymotifs.motif_atlas.chain_to_rfam_family import read_equiv_class_csv_into_dict # to get quality rank of chains
from pymotifs.motif_atlas.motifToVARNA import motif_to_varna

from fr3d.search.file_reading import readNAPairsFile
from fr3d.search.file_reading import readNAPositionsFile
from fr3d.search.search import FR3D_search, lookUpInteractions
from fr3d.search.query_processing import calculateQueryConstraints
from fr3d.search.query_processing import retrieveQueryInformation
from fr3d.search.orderBySimilarity import treePenalizedPathLength

# Modified nucleotide mappings
from fr3d.modified.mapping import modified_base_to_parent

# following variables may have to be changed depending if we are doing motif
# atlas or loop mapping

TIER1PATH          = '/usr/local/pipeline/hub-core/MotifAtlas/tier1/'
DATAPATHATLAS      = '/usr/local/pipeline/hub-core/MotifAtlas/tier1/'

DATAPATHLOOPS      = '/usr/local/pipeline/hub-core/MotifAtlas/tier2/loops'
DATAPATHRESULTS    = '/usr/local/pipeline/hub-core/MotifAtlas/tier2/search_results'
RFAMDATAPATH       = '/usr/local/pipeline/alignments'
DATAPATHALIGNMENTS = '/usr/local/pipeline/hub-core/MotifAtlas/tier2/alignments'

DATAPATHUNITS      = '/usr/local/pipeline/hub-core/data/units'  # .pickle files for each chain
DATAPATHPAIRS      = '/usr/local/pipeline/hub-core/data/pairs'  # .pickle files for each structure

SERVER = True

molecule_type = 'RNA'    # until we work on DNA

"""
All stacks and basepairs:
Note: now we have cWw and tWWa and cWB and other variations.  Make sure to deal with them right.
"""

from fr3d.classifiers.class_limits_2024 import nt_nt_cutoffs

bptypes = set()
near_bptypes = set()
for bc, cutoffs in nt_nt_cutoffs.items():
    for inter in cutoffs.keys():
        if not inter.startswith('n'):
            bptypes.add(inter)
            r = inter[0]+inter[2]+inter[1]+inter[3:]
            bptypes.add(r)
            near_bptypes.add('n'+inter)
            near_bptypes.add('n'+r)

# print('bptypes from class_limits_2024:', sorted(bptypes))
# print('near_bptypes from class_limits_2024:', sorted(near_bptypes))


# Standardizing the interactions means we don't need read class_limits_2024
bptypes = {'cWW', 'tWW', 'cWH', 'tWH', 'cWS', 'tWS', 'cHH', 'tHH', 'cHS', 'tHS',\
                   'cSS', 'tSS','cHW','tHW','cSW','tSW','cSH','tSH', 'cWB', 'cBW'}
near_bptypes = {'ntHS', 'ntSH', 'ntWH', 'ncHS', 'ncSH', 'ncWW', 'ncHH', 'ntWW', 'ntSS',\
                    'ncSS', 'ncWS', 'ntWS', 'ntHH', 'ntHW', 'ncSW', 'ncHW', 'ncWH', 'ntSW'}

stacks = {'s33','s35','s55','s53'}
near_stacks = {'ns35','ns55','ns33','ns53'}

all_stacks = stacks | near_stacks
all_bptypes = bptypes | near_bptypes

# disqualification codes for loop matches
NO_CANDIDATES          = 0
SEARCH_SPACE_CONFLICT  = -1
CONFLICTING_BASEPAIRS_AND_STACKS = 1
HAIRPIN_STACK_PENALTY1 = 3
HAIRPIN_STACK_PENALTY2 = 3
BP_PENALTY             = 4
NEAR_BP_PENALTY        = 5
STACK_PENALTY          = 6
UNMATCHED_BASEPAIR     = 7    # new 2023-06-22
FLANKING_BP_CONFLICT   = 9
MISMATCHED_BULGE       = 10

SERVER = False    # to suppress printing list lengths


def print_dictionary(d):
    for key, value in sorted(d.items()):
        print("%-20s %s" % (key, value))


# def strandify(loop_position_to_border_unit_id,all_structures=None):
#     # loops=[]
#     # for motif_group in loop_position_to_border_unit_id:
#     #     loops_alignments = motif_group["alignment"]
#     #     chainbreak = int(motif_group["chainbreak"])
#     #     for loop_id,alignment in loops_alignments.items():
#     #         pdb_id = loop_id.split("_")[1]
#     #         if pdb_id in all_structures:
#     #             new_loop = {}
#     #             new_loop["loop_id"]= loop_id
#     #             new_loop["strand"] = []
#     #             new_loop["strand"].append(alignment[0:chainbreak])
#     #             new_loop["strand"].append(alignment[chainbreak:len(alignment)])
#     #             loops.append(new_loop)
#     """
#     The code above are for json file...
#     """
#     loops = []
#     for loop_id in loop_position_to_border_unit_id:
#         fields = loop_id.split("_")
#         if fields[1] in all_structures:
#             new_loop = {}
#             new_loop["loop_id"] = loop_id
#             new_loop["strand"] = [] # will be a list of 1 strand for HL, 2 for IL, 3 for J3, etc.

#             current_strand = []
#             bordercount = 0
#             positions = sorted(loop_position_to_border_unit_id[loop_id].keys())
#             # Loop over positions in the loop and identify where strands stop and start
#             # "border" variable is 1 if the nucleotide starts or ends a single-stranded region, o/w 0
#             for position in positions:
#                 border = loop_position_to_border_unit_id[loop_id][position][0]
#                 unit_id = loop_position_to_border_unit_id[loop_id][position][1]
#                 # print("Loop %s position %s has border %s and unit id %s" % (loop_id, position, border, unit_id))
#                 current_strand.append(unit_id)
#                 bordercount += border
#                 # when you get to the second bordering nucleotide on the strand, store the strand, start a new one
#                 if bordercount == 2:
#                     bordercount = 0
#                     new_loop["strand"].append(current_strand)
#                     current_strand = []
#                 # append to the current strand until
#             loops.append(new_loop)
#     return loops


def strip_symmetry_from_pair_to_interaction(pair_to_interaction_list):
    """
    In some files, the loop unit ids do not have a symmetry operator, but the unit ids
    in pair_to_interaction has symmetry operators.
    This function strips the symmetry operator from a pair when both unit ids have the same
    symmetry operator and adds that pair to the list.
    """

    for u1,u2 in list(pair_to_interaction_list.keys()):
        f1 = u1.split("|")
        f2 = u2.split("|")

        if len(f1) == 9 and len(f2) == 9:
            # both u1 and u2 have a symmetry operator
            if f1[8] == f2[8]:
                # the symmetry operators are the same

                # remove symmetry operator and any empty fields from u1
                L = 7
                while not f1[L]:
                    L -= 1
                u1s = "|".join(f1[0:(L+1)])

                # remove symmetry operator and any empty fields from u2
                L = 7
                while not f2[L]:
                    L -= 1
                u2s = "|".join(f2[0:(L+1)])

                # store the pair with the stripped unit ids
                pair_to_interaction_list[(u1s,u2s)] = pair_to_interaction_list[(u1,u2)]

    return pair_to_interaction_list


def identify_strands(loop_position_to_border_unit_id):
    """
    Turn the data structure loop_position_to_border_unit_id into a list of dictionaries,
    each dictionary representing one loop.
    """

    loops = []
    for loop_id in loop_position_to_border_unit_id.keys():
        new_loop = {}
        new_loop["loop_id"] = loop_id
        new_loop["strand"] = [] # will be a list of 1 strand for HL, 2 for IL, 3 for J3, etc.
        new_loop["unit_ids"] = []

        if not loop_position_to_border_unit_id[loop_id]:
            # not sure why this would be None, but that can happen apparently
            continue

        if 'position' in loop_position_to_border_unit_id[loop_id]:
            # nucleotide positions have changed on this loop since the first time we did searches
            # this loop needs to be a query again and a search space again.  Rats!
            new_loop["check_again"] = True
            # delete that key
            del loop_position_to_border_unit_id[loop_id]['position']

        current_strand = []
        bordercount = 0
        positions = sorted(loop_position_to_border_unit_id[loop_id].keys())

        # Loop over positions in the loop and identify where strands stop and start
        # "border" variable is 1 if the nucleotide starts or ends a single-stranded region, o/w 0
        for position in positions:
            border = loop_position_to_border_unit_id[loop_id][position][0]
            unit_id = loop_position_to_border_unit_id[loop_id][position][1]

            # print("Loop %s position %s has border %s and unit id %s" % (loop_id, position, border, unit_id))
            current_strand.append(unit_id)
            bordercount += border

            # when you get to the second bordering nucleotide on the strand, store the strand, start a new one
            if bordercount == 2:
                bordercount = 0
                new_loop["strand"].append(current_strand)
                current_strand = []

            new_loop["unit_ids"].append(unit_id)

        loops.append(new_loop)

    return loops


def make_query_structure(loop, discrepancy = None, increasing = False):
    """
    this function should take in a loop-like-object of n strands,
    and return a query dictionary `Q` (as seen in query_definitions.py)
    Here, loop is usually the strands of the **unbulged** units of the loop
    """

    positions = get_nt_positions(loop)

    Q = defaultdict(dict) # prevents later errors on assignments

    Q['type'] = "mixed"     # geometric AND symbolic
    Q['name'] = "AVA"       # this is to DODGE Q possibly getting re-defined # could be named by loop ID

    Q['motif_atlas'] = True # when this field exists, FR3D search will not print list lengths

    if discrepancy:
        Q['discrepancy'] = discrepancy
    else:
        if len(loop['strand']) <= 2:
            Q['discrepancy'] = 1.0  # 1.0 is what the Motif Atlas uses for HL and IL
        else:
            Q['discrepancy'] = 1.0  # in case we want to try something new for J3

    Q = make_num_positions(Q, positions)

    Q["interactionMatrix"] = emptyInteractionMatrix(Q['numPositions'])
    Q['unitID'] = []        # these are the unit ids that the query will actually use for the search
    Q['requiredUnitType'] = [None] * Q['numPositions']
    Q['requiredMoleculeType'] = [molecule_type] * Q['numPositions']

    Q['errorMessage'] = []
    Q['userMessage'] = []
    Q['MAXTIME'] = 60*60    # maximum search time 1 hour, may not be imposed
    Q["CPUTimeUsed"] = 0    # needs to be tracked, but actually ignored here

    Q = add_FR3D_diagnostics(Q)

    # below are things that prefer to be fed strand by strand
    for index, strand in enumerate(loop['strand']):
        Q = make_unit_ids(Q, strand)
        if len(positions[index]) > 1:
            # if strand is more than 1 nt, find direction + constraint, else `pass`
            Q = make_directional_constraints(Q, positions[index], increasing)

    # below are things requiring above functions to be done, but dislike strand by strand
    Q = make_search_files(Q)

    if len(positions) == 0:
        print("make_query_structure: No nt positions in loop %s" % loop)
        Q['errorStatus'] == "write and exit"
        return Q
    elif len(positions[0]) == 0:
        print("make_query_structure: No nt positions in first strand of loop %s" % loop)
        Q['errorStatus'] == "write and exit"
        return Q

    Q = make_req_interactions(Q, positions)

    return Q


def add_FR3D_diagnostics(Q):
    # turn on lots of FR3D output for diagnostics

    # Q["printInputQuery"] = 0
    # Q["printProcessedQuery"] = 0
    # Q["printSearchFiles"] = 10
    # Q["printSearchingFile"] = 1
    # Q["printNotEnoughUnits"] = 1
    # Q["printListLengths"] = 1
    # Q["printFoundCandidates"] = 1
    # Q["printFoundPossibilities"] = 1
    # Q["printTimer"] = 1

    return Q

def make_flanking_bp_query_structure(loop):
    """
    A small query for IL and J3 and J4 flanking basepairs; not for hairpins.
    this function should take in a loop-like-object of n strands,
    and return a query dictionary `Q` (as seen in query_definitions.py)
    """

    positions = get_nt_positions(loop)
    flanking_positions = [p for p in positions]

    Q = defaultdict(dict) # prevents later errors on assignments

    Q["DATAPATHUNITS"] = DATAPATHUNITS
    Q["DATAPATHPAIRS"] = DATAPATHPAIRS
    Q["activeInteractions"] = ["cWW"]
    Q["PDB_data_file"] = []    # block fr3d.search code from adding it
    Q["MAXTIME"] = 60*60
    Q["CPUTimeUsed"] = 0

    Q = add_FR3D_diagnostics(Q)

    if len(positions) == 0:
        print("No nt positions in loop %s" % loop)
        Q['errorStatus'] == "write and exit"
        return Q
    elif len(positions[0]) == 0:
        print("No nt positions in first strand of loop %s" % loop)
        Q['errorStatus'] == "write and exit"
        return Q

    for i in range(len(positions)):
        flanking_positions[i] = [ positions[i][0],positions[i][-1]]

    # flanking nt positions
    flanking_strands = [s for s in loop["strand"]]
    for i in range(len(loop['strand'])):
        flanking_strands[i] = [ loop['strand'][i][0],loop['strand'][i][-1]]

    # philosophy = if the key,value takes more than one line to assign, define a `make` function for it
    Q['type'] = "mixed"     # geometric AND symbolic
    Q['name'] = "AVA"       # this is to DODGE Q possibly getting re-defined # could be named by loop ID
    Q['motif_atlas'] = True # make it possible to tailor output for motif atlas AVA

    if len(loop['strand']) <= 2:
        Q['discrepancy'] = 1.0  # 1.0 is what the Motif Atlas uses for HL and IL
    else:
        Q['discrepancy'] = 1.0  # in case we want to try something new for J3

    Q = make_num_positions(Q, flanking_positions)

    Q["interactionMatrix"] = emptyInteractionMatrix(Q['numPositions'])

    # add ncWW below because some old loops are now annotated with ncWW flanking pairs

    if Q['numPositions'] == 4:
        # Only for IL
        Q["interactionMatrix"][0][3] = "cWW ncWW and GC CG AU UA GU UG +mod"
        Q["interactionMatrix"][1][2] = "cWW ncWW and GC CG AU UA GU UG +mod"
        Q["interactionMatrix"][1][0] = ">"
        Q["interactionMatrix"][3][2] = ">"
    elif Q['numPositions'] > 4 and Q['numPositions'] % 2 == 0:
        # J3, J4, etc.
        N = Q['numPositions']

        Q["interactionMatrix"][0][N-1] = "cWW ncWW and GC CG AU UA GU UG +mod"
        for i in range(1,N-2,2):
            Q["interactionMatrix"][i][i+1] = "cWW ncWW and GC CG AU UA GU UG +mod"

        for i in range(1,N,2):
            Q["interactionMatrix"][i][i-1] = ">"

        # Q["interactionMatrix"][1][2] = "cWW and GC CG AU UA GU UG"
        # Q["interactionMatrix"][3][4] = "cWW and GC CG AU UA GU UG"
        # Q["interactionMatrix"][1][0] = ">"
        # Q["interactionMatrix"][3][2] = ">"
        # Q["interactionMatrix"][5][4] = ">"

    Q['unitID'] = []
    Q['requiredUnitType'] = [None] * Q['numPositions'] # these list operations are copied from query_processing.py
    Q['requiredMoleculeType'] = [molecule_type] * Q['numPositions']

    Q['errorMessage'] = []
    Q['userMessage'] = []

    # below are things that prefer to be fed strand by strand
    for index, strand in enumerate(flanking_strands):
        Q = make_unit_ids(Q, strand)

    # below are things requiring above functions to be done, but dislike strand by strand
    Q = make_search_files(Q)

    return Q


def get_nt_positions(loop):
    """
    takes in a list of strings representing nucleotides
    returns a list of lists of nucleotide position numbers (by strand)
    ex: positions = [[nt1, nt2, nt3], [nt4, nt5, nt6]]
    ex: positions = [[0, 1, 2], [3, 4, 5]]
    ^ nts 1-3 are on one strand and nt4-6 are on another
    also, nts here look like: "5J7L|1|AA|C|569"
    """

    positions = []
    count = 0
    for index, strand in enumerate(loop['strand']):
        tempList = []
        for element in strand:
            tempList.append(count)
            count += 1
        positions.append(tempList)

    return positions


def make_directional_constraints(query, positions, increasing = False):
    """
    Create FR3D query constraints to indicate that nucleotides
    on each strand must be in increasing chain index order.
    That focuses the FR3D search and speeds it up.
    """

    if increasing:
        # all positions must be in 5' to 3' order, even across strand breaks
        # not good for IL, which can match in either rotation
        # reasonable for J3 and higher, where the combinatorics of many strands
        # make it nearly impossible to match in a different rotation
        for i in range(1,query['numPositions']):
            query['interactionMatrix'][i][i-1] = ">"
    else:
        # all positions within each strand must be in 5' to 3' order
        previousPositions = [] # needed to FILL lower diagonal
        firstIter = True
        for position in positions: # remember that this function is called on each strand
            if firstIter == False: # if first iteration, previousPositions is empty
                for prev in previousPositions:
                    query['interactionMatrix'][position][prev] = ">"
                previousPositions.append(position)
            else:
                firstIter = False
                previousPositions.append(position)

    return query


def make_num_positions(query, positions):
    """
    counts the nts in a loop
    """

    counter = 0
    for strand in positions:
        counter += len(strand)
    query['numPositions'] = counter

    return query


def emptyInteractionMatrix(n):
    """
    set up an empty list of required interactions,
    then the user only needs to set the actual constraints needed
    """

    emptyList = defaultdict(dict)
    for i in range(0, n):
        for j in range(0, n):
            emptyList[i][j] = ""

    return emptyList


def make_req_interactions(query, positionsList):
    """
    fills out a matrix of required cWWs:
    last nt of first strand needs to cWW to first nt of second strand (and so on)
    should be putting "cWW" ABOVE diagonal
    """

    heads = [] #remember that lists are ordered, so these are parallel
    tails = []
    for strand in positionsList: # filling points of interest
        heads.append(strand[0])
        tails.append(strand[-1])

    prevTail = -1
    firstHead = -1
    for head, tail in zip(heads, tails): # adding constraint
        # first tail -> second head, second tail -> third head,
        # ... last tail -> first head
        if firstHead == -1: # if first pass
            firstHead = head
        else: # if not first pass, we have connections to make
            query['interactionMatrix'][prevTail][head] = "cWW ncWW GC CG AU UA GU UG +mod"
        if tail == tails[-1]: # if last pass
            query['interactionMatrix'][firstHead][tail] = "cWW ncWW GC CG AU UA GU UG +mod"
        prevTail = tail

    return query


def make_search_files(query):
    """
    this function assumes `make_unit_ids()` was already run
    Q['searchFiles'] = [FIRST part of unit_id]
    """

    chains = set()
    for nt in query['unitID']:
        fields = nt.split("|")
        link = "|".join(fields[0:3])
        chains.add(link)
    query['searchFiles'] = ["+".join(chains)] # future, don't overwrite

    return query


def make_unit_ids(query, strand):
    """
    example of what this should make
    Q["unitID"] = ["4V9F|1|0|U|1026","4V9F|1|0|A|1032","4V9F|1|0|G|1034", ... ]
    """

    for nt in strand:
        query['unitID'].append(nt)

    return query


def combine_dicts(x, y):
    z = x.copy()
    z.update(y)
    return(z)


def readPAI(query, ifename, alternate = "", by_chain = False):
    """
    WAS readPositionsAndInteractions() from ifedata.py.
    It got copied + changed here for AvA needs
    """

    NA_positions_file_name = ifename.replace("|","-").replace("+","_")
    starting_index = 0
    ifedata = {}
    ifedata['index_to_id'] = {}
    ifedata['id_to_index'] = {}
    ifedata['centers'] = np.empty((0, 3))
    ifedata['rotations'] = np.empty((0, 3, 3))
    ifedata['ids'] = []
    ifedata['units'] = []
    allModels = []

    for chain_string in NA_positions_file_name.split("_"):
        # print("in readPAI reading %s" % (chain_string))
        query, centers, rotations, ids, id_to_index, index_to_id, chainIndices = readNAPositionsFile(query, chain_string, starting_index)

        if len(centers) > 0:
            # index is the internal number of the nucleotide in FR3D, not the chain index
            ifedata['index_to_id'] = combine_dicts(ifedata['index_to_id'], index_to_id)
            ifedata['id_to_index'] = combine_dicts(ifedata['id_to_index'], id_to_index)
            ifedata['centers'] = np.append(ifedata['centers'], centers, axis = 0)
            ifedata['rotations'] = np.append(ifedata['rotations'], rotations, axis = 0)
            ifedata['ids'].extend(ids)

            for unitID in ids:
                unit_information = {}
                unit_information["centers"] = centers[id_to_index[unitID]-starting_index]
                unit_information["rotations"] = rotations[id_to_index[unitID]-starting_index]
                data = unitID.split("|")
                allModels.append(data[1])                 # extract model number
                unit_information["unitType"] = data[3]    # get the nucleotide sequence
                unit_information["moleculeType"] = molecule_type  # could be DNA eventually ...
                unit_information["chainindex"] = chainIndices[id_to_index[unitID]-starting_index]

                # allow flanking cWW pair constraints to match to modified bases
                p = modified_base_to_parent.get(data[3], data[3])
                unit_information["parentType"] = p
                if p in ['DA','DC','DG']:
                    p = p[1]
                elif p == 'DT':
                    p = 'U'
                unit_information["parentAsRNA"] = p

                ifedata["units"].append(unit_information)

            starting_index += len(centers)

    PDBID = ifename.split("|")[0]
    query, interactionToPairs, pairToInteractions, pairToCrossingNumber = readNAPairsFile(query, PDBID, ifedata["id_to_index"], alternate)
    ifedata['interactionToPairs'] = interactionToPairs
    ifedata['pairToInteractions'] = pairToInteractions
    ifedata['pairToCrossingNumber'] = pairToCrossingNumber
    ifedata['models'] = allModels

    return query, ifedata


def create_query(loop, ifedata = {}, discrepancy = None, increasing = False):
    """
    Make a query data structure based on the data in loop
    Makes the query (Q) without bulges, while bulges should be left in search spaces (ifedata)
    """

    loop_type = loop["loop_id"].split("_")[0]

    if not ifedata:
        print("No ifedata in %s" % loop["loop_id"])
        return None

    query_length = sum([len(strand) for strand in loop["strand"]])
    query_bulged_length = len(loop["bulged"])

    # identify non-bulged nucleotides to form the basis of the query
    # double list comprehension to keep the structure of ['strand'] being a list of lists, while removing bulges
    unbulged = {}
    unbulged['strand'] = [[unitID for unitID in strand if unitID not in loop['bulged']] for strand in loop['strand']]

    added_nts = False
    if loop_type == "HL" and len(unbulged['strand']) > 0 and len(unbulged['strand'][0]) == 2:
        # add nucleotides to search with
        print('%s needs more nucleotides, has unbulged: %s' % (loop["loop_id"],unbulged['strand']))
        strand = loop['strand'][0]
        if len(strand) >= 3:
            unbulged['strand'] = [[strand[0], strand[1], strand[-1]]]
            print('%s getting more nucleotides: %s' % (loop["loop_id"],unbulged['strand']))
            added_nts = True
            if strand[1] in loop["bulged"]:
                loop["bulged"].remove(strand[1])

    # when a 5-nt IL with one bulged nt, keep all 5 nucleotides in unbulged variable
    if loop_type == "IL" and query_length == 5 and query_bulged_length == 1:
        query = make_query_structure(loop,98)
    else:
        query = make_query_structure(unbulged,discrepancy,increasing)

    if 'errorStatus' in query and query['errorStatus'] == "write and exit":
        return None

    query['fullUnits'] = [unitID for strand in loop['strand'] for unitID in strand]

    # if we already have centers and rotations
    if query_bulged_length == 0 or (query_length == 5 and query_bulged_length == 1):
        query['centers'] = ifedata['centers']
        query['rotations'] = ifedata['rotations']
    else:
        # do not store centers or rotations for bulged nucleotides
        indices_of_unbulged = []
        for unit_id in query['fullUnits']:
            if not unit_id in loop['bulged']:
                indices_of_unbulged.append(ifedata['id_to_index'][unit_id])

        query['centers'] = ifedata['centers'][indices_of_unbulged]
        query['rotations'] = ifedata['rotations'][indices_of_unbulged]

    if 'errorStatus' in query and query['errorStatus'] == "write and exit":
        return None
    else:
        query = calculateQueryConstraints(query)
        query['numStrands'] = len(loop['strand'])

    query['added_nts'] = added_nts

    # Starting on 2/6/2025, explicitly ask to tolerate "flips" in orientation error
    if loop_type == 'IL':
        strand0unbulged = [u for u in loop['strand'][0] if not u in loop['bulged']]
        strand1unbulged = [u for u in loop['strand'][1] if not u in loop['bulged']]

        if len(strand0unbulged) == 3 and len(strand1unbulged) == 3:
            # except don't do that when the query is a 3x3 IL
            query['3x3'] = True
            print(strand0unbulged,' 3x3 ',strand1unbulged)
            query["flip"] = False   # added for 3.360
        else:
            query["flip"] = True
    else:
        query["flip"] = True

    return query


def create_flanking_bp_query(loop, ifedata = {}):

    unbulged = {}

    unbulged['strand'] = [[unitID for unitID in strand if unitID not in loop['bulged']] for strand in loop['strand']]
    flanking_bp_query = make_flanking_bp_query_structure(unbulged)

    if not flanking_bp_query:
        return None

    flanking_bp_query = retrieveQueryInformation(flanking_bp_query)

    if 'errorStatus' in flanking_bp_query and flanking_bp_query['errorStatus'] == "write and exit":
        return None
    elif not flanking_bp_query:
        return None
    else:
        flanking_bp_query = calculateQueryConstraints(flanking_bp_query)
        flanking_bp_query['numStrands'] = len(loop['strand'])

    return flanking_bp_query


def get_bulge_old(query,search_space):
    query_bulge = ''
    search_space_bulge = ''
    if len(query['loop_info']['bulged'])>0:
        query_bulge = query['loop_info']['bulged'][0].split("|")[3]
    if len(search_space['loop_info']['bulged'])>0:
        search_space_bulge = search_space['loop_info']['bulged'][0].split("|")[3]
    return(query_bulge,search_space_bulge)


def get_bulge(query):
    """
    For single bulge IL, return the parent of the bulged base
    Also return the position of the bulge in the loop
    """

    bulge = ""
    if len(query['loop_info']['bulged']) > 0:
        bulge_unit_id = query['loop_info']['bulged'][0]
        bulge = bulge_unit_id.split("|")[3]

    parent = modified_base_to_parent.get(bulge, bulge)

    return parent


def filter_on_unmatched_nucleotides(candidates, search_space, query):
    """
    The FR3D search finds nucleotides in the search space that match the nucleotides in the query.
    But some nucleotides in the search space might not be matched, maybe they are bulged out, maybe another reason.
    Analyze unmatched nucleotides of candidates that have no disqualification yet
    Check interactions between unmatched nucleotides first,
    Then interactions between unmatched and matched nucleotides
    """

    candidates_data = search_space['ifedata']

    not_rejected = []

    for candidate in candidates:
        if len(candidate["dq"]) > 0:
            # already will be disqualified, no need to check further
            continue

        matched_indices = candidate["indices"]
        unmatched_nts = list(set(candidates_data["index_to_id"].keys()) - set(matched_indices))
        if len(unmatched_nts) == 0:
            not_rejected.append(candidate)
            continue

        stack_counter = 0
        for i in range(len(unmatched_nts)):
            # Look for interactions between unmatched nucleotides
            for j in range(i+1,len(unmatched_nts)):
                interaction_with_unmatched = get_set_of_interactions(unmatched_nts[i],unmatched_nts[j],candidates_data)
                if interaction_with_unmatched.intersection(bptypes):
                    candidate['dq'].append(BP_PENALTY)
                if interaction_with_unmatched.intersection(near_bptypes) :
                    candidate['dq'].append(NEAR_BP_PENALTY)

            # Look for interaction between matched and unmatched nucleotides
            for index in matched_indices:
                interaction_with_matched = get_set_of_interactions(unmatched_nts[i],index,candidates_data)
                if interaction_with_matched.intersection(bptypes) :
                    candidate['dq'].append(BP_PENALTY)
                if interaction_with_matched.intersection(near_bptypes):
                    candidate['dq'].append(NEAR_BP_PENALTY)
                if interaction_with_matched.intersection(stacks):
                    stack_counter += 1

        if stack_counter > 1:
            candidate['dq'].append(STACK_PENALTY)

        candidate['dq'] = list(set(candidate['dq']))

        if len(candidate['dq']) == 0:
            not_rejected.append(candidate)

    return not_rejected


def filter_out_conflicting_basepairs_and_stacks(candidates, query, search_space, query_unit_id=None):
    """
    Look for conflicting base pairs or stacking between query and candidates
    Return candidates with no conflicts
    candidates data is just the ifedata of the candidate
    Data needed from query: ifedata, Q['unitID']
    Some loop mapping queries use bulged bases and so cannot rely on query["Q"]["unitID"]
    """

    query_ifedata = query['ifedata']
    if not query_unit_id:
        query_unit_id = query["Q"]["unitID"]

    candidates_data = search_space['ifedata']

    not_rejected = []

    # sort candidates by discrepancy
    candidates = sorted(candidates, key = lambda i: i['discrepancy'])

    for candidate in candidates:

        # record disqualification code(s)
        DQ_code = []

        filtered_candidate = {}
        matched_search_space_indices = set()

        # loop over pairs of aligned positions, compare the basepairs and stacks there
        for i in range(len(candidate['indices'])): # goes like 0 to 7
            # next line utilizes Q not including bulges, and how its 'unitID's are in index order, to
            # return to its 'ifedata's indexing system
            index_i_of_query = query_ifedata['id_to_index'][query_unit_id[i]]
            index_i_of_candidate = candidate['indices'][i] #The integer indentifying nucleotide in i position
            matched_search_space_indices.add(index_i_of_candidate)

            for j in range(i + 1, len(candidate['indices'])): # goes like i+1 to 7
                index_j_of_query = query_ifedata['id_to_index'][query_unit_id[j]]
                index_j_of_candidate = candidate['indices'][j]

                query_interaction = get_set_of_interactions(index_i_of_query, index_j_of_query, query_ifedata) & (bptypes | stacks)
                candidate_interaction = get_set_of_interactions(index_i_of_candidate, index_j_of_candidate, candidates_data) & (bptypes | stacks)
                DQ_code += are_conflicting(query_interaction, candidate_interaction)

        if len(DQ_code) == 0:
            filtered_candidate['dq'] = []
            filtered_candidate['query_id'] = query['loop_info']['loop_id']
            filtered_candidate['search_space_id'] = search_space['loop_info']['loop_id']
            filtered_candidate['indices'] = candidate['indices']
            filtered_candidate['query_unit_ids'] = query_unit_id
            filtered_candidate['target_unit_ids'] = candidate['unitids']
            filtered_candidate['discrepancy'] = candidate['discrepancy']
            filtered_candidate['numStrands'] = query['Q']['numStrands']

            # filtered_candidate['query_rfam'] = query['loop_info']['rfam']
            # filtered_candidate['search_space_rfam'] = search_space['loop_info']['rfam']

        else: # if has a real DQ_code
            filtered_candidate['dq'] = list(set(DQ_code))
            filtered_candidate['discrepancy'] = candidate['discrepancy']
            filtered_candidate['numStrands'] = query['Q']['numStrands']

        not_rejected.append(filtered_candidate)

    return not_rejected


def keep_lowest_discrepancy_candidate(candidates):
    """
    this takes a list of possible candidates and returns the one
    with the lowest discrepancy, unless all candidates have
    disqualification codes, then it returns an empty list
    """
    lowest_dq_discrepancy = 100 # purposely seeded high
    lowest_matched_discrepancy = 100
    best_match = 0
    best_dq = None

    for candidate in candidates:
        if len(candidate['dq']) == 0:
            if candidate['discrepancy'] < lowest_matched_discrepancy:
                lowest_matched_discrepancy = candidate['discrepancy']
                best_match = candidate
        else: # if a disqualified match
            if candidate['discrepancy'] < lowest_dq_discrepancy:
                lowest_dq_discrepancy = candidate['discrepancy']
                best_dq = candidate

    if best_match:
        return [best_match]
    elif best_dq:
        # if all disqualified, return the best so we can note the disqualification
        return [best_dq]
    else:
        return []


def get_set_of_interactions(i, j, ifedata):
    """
    the definition of interaction here is hard stacks or hard basepairs,
    but not near stacks nor near basepairs
    """
    if (i,j) in ifedata['pairToInteractions']:
        return set(ifedata['pairToInteractions'][(i,j)])
    else:
        return set([])


def are_conflicting(q_int, c_int):
    """
    q_int is short for "interaction of query between nts (i, j)"
    c_int is short for "interaction of canditate between nts (i, j)"
    the definition of interaction here is hard stacks or hard basepairs,
    but not near stacks nor near basepairs
    """
    # return options [-1 = "compatible", 8 = "basestack_mismatch",
    #                   7 = "basepair_mismatch"]
    q_int = set(q_int)
    c_int = set(c_int)
    BASESTACK_MISMATCH = 8
    BASEPAIR_MISMATCH = 7

    if q_int == c_int:
        return []
    if (q_int & stacks) and (c_int & stacks): # loose due to modeling flips
        return []
    if not ( q_int&stacks) and (c_int&stacks):
        if len(q_int) != 0:
            return [BASESTACK_MISMATCH] # base stack mismatch
        else:
            return []
    if (q_int & stacks) and not (c_int & stacks):
        if len(c_int) != 0:
            return [BASESTACK_MISMATCH] # base stack mismatch
        else:
            return []
    # no stacks by this line
    if len(q_int) == 0 or len(c_int) == 0: # one is empty
        return [] # compatible
    if q_int != c_int: # i dont believe any loops have made it this far yet (82x82 tested)
        return [BASEPAIR_MISMATCH] # base pair mismatch
    print("we hit the default case... inspect c_int and q_int")
    raise BaseException("error in `are_conflicting()` logic.\nq_int = " + q_int + "\nc_int = " + c_int)


def name_structure(dict_of_searches):
    """
    this is a helper function for the naming of files
    """
    number_of_strands = dict_of_searches[list(dict_of_searches.keys())[0]][0]['numStrands']

    if number_of_strands == 1:
        return("HL") # hairpin loop
    if number_of_strands == 2:
        return("IL") # internal loop
    if number_of_strands > 2:
        return("J" + str(number_of_strands))
    else:
        raise BaseException("number of strands is {}".format(number_of_strands))


def load_previous_search_results_one_pdb(loop_type,pdb_id, path = DATAPATHRESULTS):

    results = {}
    if not os.path.exists(path):
        os.mkdir(path)
        return results

    file_name = loop_type + "_" + pdb_id + "_search_results.pickle.gz"
    file_path = os.path.join(path,file_name)

    if os.path.exists(file_path):
        with gzip.open(file_path, 'rb') as f:
            results = pickle.load(f)

        for (query,search_space),value in results.items():
            # Load No-match result as [{'dq': 0, 'discrepancy': 99}] instead of an integer

            # if query == 'J3_5J7L_036':
            #     print(query,search_space,value)

            if isinstance(value,int):
                results[(query,search_space)] = [{'dq': [value], 'discrepancy': 99}]

    return results


def save_search_results_one_pdb(loop_type, pdb_id, results, path = DATAPATHRESULTS):
    if not os.path.exists(path):
        os.mkdir(path)

    file_name = loop_type + '_' + pdb_id + '_search_results.pickle.gz'
    file_path = os.path.join(path, file_name)

    if sys.version_info[0] < 3:
        with gzip.open(file_path,'wb') as f:
            pickle.dump(results, f)
    else:
        with gzip.open(file_path,'wb') as f:
            pickle.dump(results, f, protocol = 4)

    return


def list_all_chains_in(loops_from_structure, best_quality_chains = False):
    from itertools import permutations
    accumulator_of_chains = []

    # this top "dict" block... I'm not sure if it's used anymore (09/18/2023 - Adam)
    if type(loops_from_structure) is dict:
        for loop_id, loop_info in loops_from_structure.items():
            this_loop_chains = []
            this_loop_unit_ids = [nt for strand in loop_info['loop_info']['strand'] for nt in strand]
            for unit_id in this_loop_unit_ids:
                fields = unit_id.split("|")
                one_chain = "|".join(fields[0:3])
                this_loop_chains.append(one_chain)
            accumulator_of_chains = accumulator_of_chains + ["+".join(set(sorted(this_loop_chains)))]

    # if "loop_info" in loops_from_structure:
    if type(loops_from_structure) is list:
        for loop_dict in loops_from_structure:
            this_loop_chains = []
            this_loop_unit_ids = [nt for strand in loop_dict['strand'] for nt in strand]
            for unit_id in this_loop_unit_ids:
                fields = unit_id.split("|")
                one_chain = "|".join(fields[0:3])
                this_loop_chains.append(one_chain) # HOW CAN WE KEEP IFEs IN ORDER?
            if len(set(this_loop_chains)) == 1: # begin new lines
                accumulator_of_chains.append(this_loop_chains[0])
            else: # append all combinations  "chain1+chain2" AND "chain2+chain1"
                for ife_possibility in permutations(set(this_loop_chains)):
                    accumulator_of_chains.append("+".join(ife_possibility))

            # old method instead of the above if/else
            # accumulator_of_chains = accumulator_of_chains + ["+".join(set(this_loop_chains))]

    accumulator_of_chains = list(set(accumulator_of_chains))

    if best_quality_chains:
        best_chains = {} # eventually: best_chains[equiv_class] = {chain = 'XYZ|1|1', quality_score = '1'}
        chain_to_equiv_class = read_equiv_class_csv_into_dict() # secondary keys = 'equivalence class', 'quality rank'
        for chain in accumulator_of_chains:
            # print(f"on chain {chain}... equiv: {chain_to_equiv_class[chain]['equivalence class']}... quality_rank: {chain_to_equiv_class[chain]['quality rank']}")
            if chain_to_equiv_class.get(chain): # if key exists
                equiv_class = chain_to_equiv_class[chain]['equivalence class']
                quality_rank = chain_to_equiv_class[chain]['quality rank']
                if not best_chains.get(equiv_class):
                    best_chains[equiv_class] = {'chain': chain, 'quality_score': quality_rank}
                else: # if we have another chain from this equiv class
                    if quality_rank < best_chains[equiv_class]['quality_score']:
                        # we found a better chain. Overwrite
                        best_chains[equiv_class] = {'chain': chain, 'quality_score': quality_rank}
            # else:
            #     # keep? yes for now
            #     best_chains[chain] = {'chain': chain}

        # now we have the best chains in a nested dict. Time to retrieve them
        accumulator_of_chains = [] # reset
        for equiv_class, info in best_chains.items():
            accumulator_of_chains.append(info['chain'])

    return(accumulator_of_chains)


def nt_to_loop_mapping(search_spaces):
    """
    this is a needed dict for later trying to map nucleotides of the search space to the
    query's loop to search it with
    """

    mapping = {}

    for loop_name_in_ss, search_space_loop_data in search_spaces.items():
        # this next line is a list comprehension, to flatten a list of lists into a regular 1D list
        nts_of_ss_loop = [nt for strand in search_space_loop_data['loop_info']['strand'] for nt in strand]

        for nt in nts_of_ss_loop:
            mapping[nt] = loop_name_in_ss


    return(mapping)


def add_bulged_nucleotides_to_loop(loop, inter_lst):
    """
    Analyze one loop and identify the bulged nucleotides
    {id: HL_4V9F_003, strands: [[4V9F|1|0|G|421,4V9F|1|0|G|422,4V9F|1|0|C|423]],bulged: [4V9F|1|0|G|422]}
    """
    non_bulged_interactions = all_bptypes | all_stacks

    bulged_nts = []
    all_nts = [unit_id for sublist in loop['strand'] for unit_id in sublist]

    # in some cases, the flanking pair is very bad and flanking nucleotides
    # are not making an interaction, and so were being marked as bulged!
    interior_unit_ids = []
    # loop over strands
    for strand in loop['strand']:
        # loop over second to second to last nucleotide in the strand
        for unit1 in strand[1:-1]:
            interior_unit_ids.append(unit1)

    # print(loop,interior_unit_ids,all_nts)

    for unit1 in interior_unit_ids:

        unit1_interactions = set()

        for unit2 in all_nts:
            unit1_interactions.update(inter_lst[(unit1, unit2)])

        if not unit1_interactions.intersection(non_bulged_interactions):
            bulged_nts.append(unit1)

    return bulged_nts


def create_queries_and_search_spaces(loops,discrepancy=None):

    queries = {}
    flanking_bp_queries = {}
    search_spaces = {}
    pdb_data = {}
    current_pdb = ''

    # old_total_memory = 0

    # make sure loops is sorted by PDB id, for efficient reading of interactions
    for loop in sorted(loops, key=lambda k: k['loop_id']):
        # get a flat list of the nucleotides (full unit ids) of this loop
        this_loop_unit_ids = [nt for strand in loop['strand'] for nt in strand]

        # print('Creating search space and query for loop %s with %2d nucleotides' % (loop['loop_id'], len(this_loop_unit_ids)))

        # print total memory usage
        # total_memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        # print('Memory usage: %s (kb) delta %s' % (total_memory, total_memory - old_total_memory))
        # old_total_memory = total_memory

        # print memory used by variable queries
        # print('Memory usage of queries: %s (kb)' % sys.getsizeof(queries))

        # print all active variables and the memory they use; I don't know which variable is using all the memory
        # for var, obj in sorted(locals().items()):
        #     s = sys.getsizeof(obj)
        #     if s > 500:
        #         print('%8d %s' % (s,var))

        # print all global variables and the memory they use
        # for var, obj in sorted(globals().items()):
        #     s = sys.getsizeof(obj)
        #     if s > 500:
        #         print('%8d %s' % (s,var))

        pdb = loop['loop_id'].split("_")[1]

        # if we're on a new pdb, read all pairs
        if pdb != current_pdb:
            pdb_data = {}
            current_pdb = pdb

            # read all pairwise interactions in this pdb structure
            pair_to_interaction_list, pair_to_crossing_number = get_fr3d_pair_to_interaction_list(pdb, near = True, standard = True)

            pair_to_interaction_list = strip_symmetry_from_pair_to_interaction(pair_to_interaction_list)

        # add bulged nucleotide data to this loop
        loop["bulged"] = add_bulged_nucleotides_to_loop(loop, pair_to_interaction_list)

        this_loop_chains = set()
        for unit_id in this_loop_unit_ids:
            fields = unit_id.split("|")
            one_chain = "|".join(fields[0:3])
            this_loop_chains.add(one_chain)

        for chain in this_loop_chains:
            if not chain in pdb_data:
                Q = {}
                Q['DATAPATHUNITS'] = DATAPATHUNITS
                Q['DATAPATHPAIRS'] = DATAPATHPAIRS
                Q['activeInteractions'] = all_bptypes | all_stacks
                Q['errorMessage'] = []
                Q['userMessage'] = []
                Q, pdb_data[chain] = readPAI(Q, chain, by_chain = True)

        ifedata = {}
        ifedata['index_to_id'] = {}
        ifedata['id_to_index'] = {}
        ifedata['centers'] = np.empty((0, 3))
        ifedata['rotations'] = np.empty((0, 3, 3))    # is this actually needed?
        ifedata['ids'] = this_loop_unit_ids
        ifedata['units'] = []                         # heavily used by FR3D search
        ifedata['interactionToPairs'] = {}            # key will be interaction type
        ifedata['pairToInteractions'] = defaultdict(list)
        ifedata['pairToCrossingNumber'] = {}
        ifedata['models'] = []

        # convert chain_indices from the full pdb_data
        # to shallow_indices (loop indices), and extract only data
        # relevant to this loop
        for shallow_index, unit_id in enumerate(this_loop_unit_ids):
            fields = unit_id.split("|")
            chain = "|".join(fields[0:3])

            data_index = None  # find where this unit_id is in the pdb_data

            if unit_id in pdb_data[chain]['id_to_index']:
                data_index = pdb_data[chain]['id_to_index'][unit_id]
            else:
                # Sometimes unit id like 1NUV|1|A|G|26 needs to match 1NUV|1|A|G|26||A,
                # or the other way around. Try both.
                fields = unit_id.split("|")
                if len(fields) == 5:
                    try_id = unit_id + "||A"
                elif len(fields) == 7 and fields[6] == "A":
                    try_id = "|".join(fields[0:5])
                else:
                    try_id = None
                if try_id and try_id in pdb_data[chain]['id_to_index']:
                    data_index = pdb_data[chain]['id_to_index'][try_id]

            if data_index == None:
                print("create_queries_and_search_spaces: Could not find internal data index for %s in %s ====" % (unit_id,loop['loop_id']))
                ifedata['missing_center'] = True
            else:
                ifedata['id_to_index'][unit_id] = shallow_index
                ifedata['index_to_id'][shallow_index] = unit_id
                ifedata['centers'] = np.append(ifedata['centers'], np.asarray([pdb_data[chain]['centers'][data_index]]), axis = 0)
                ifedata['rotations'] = np.append(ifedata['rotations'], np.asarray([pdb_data[chain]['rotations'][data_index]]), axis = 0)
                ifedata['units'].append(pdb_data[chain]['units'][data_index])
                ifedata['models'].append(pdb_data[chain]['models'][data_index])

        # get the required interactions from interactionToTriple
        for shallow_index_1, unit_id_1 in enumerate(this_loop_unit_ids):
            for shallow_index_2, unit_id_2 in enumerate(this_loop_unit_ids):
                if (unit_id_1, unit_id_2) in pair_to_interaction_list:
                    # store the interaction according to the shallow index, not just the unit ids
                    index_pair = (shallow_index_1, shallow_index_2)

                    ifedata['pairToInteractions'][index_pair] = pair_to_interaction_list[(unit_id_1,unit_id_2)]
                    ifedata['pairToCrossingNumber'][index_pair] = pair_to_crossing_number[(unit_id_1,unit_id_2)]

                    crossing_number = pair_to_crossing_number[(unit_id_1,unit_id_2)]

                    # needed by FR3D search
                    for interaction in pair_to_interaction_list[(unit_id_1, unit_id_2)]:
                        if interaction in all_bptypes | all_stacks:
                            if not interaction in ifedata['interactionToPairs']:
                                ifedata['interactionToPairs'][interaction] = [[],[]]
                            ifedata['interactionToPairs'][interaction][0].append(index_pair)
                            ifedata['interactionToPairs'][interaction][1].append(crossing_number)

        if not 'missing_center' in ifedata:

            Q = create_query(loop, ifedata, discrepancy)

            if not Q:
                continue
            elif "errorStatus" in Q and Q["errorStatus"] == "write and exit":
                continue

            if not loop["loop_id"].startswith('HL'):
                bp_Q = create_flanking_bp_query(loop, ifedata)

                # size of bp_Q
                # print('bp_Q',sys.getsizeof(bp_Q))

                if not bp_Q:
                    continue

                if 'PDB_data_file' in bp_Q and len(bp_Q['PDB_data_file']) > 0:
                    del bp_Q['PDB_data_file']

                flanking_bp_queries[loop['loop_id']] = {}
                flanking_bp_queries[loop['loop_id']]['loop_info'] = loop
                flanking_bp_queries[loop['loop_id']]['Q'] = bp_Q
                flanking_bp_queries[loop['loop_id']]['ifedata'] = ifedata

            queries[loop['loop_id']] = {}
            queries[loop['loop_id']]['loop_info'] = loop
            queries[loop['loop_id']]['Q'] = Q
            queries[loop['loop_id']]['ifedata'] = ifedata

            search_spaces[loop['loop_id']] = {}
            search_spaces[loop['loop_id']]['ifedata'] = ifedata
            search_spaces[loop['loop_id']]['loop_info'] = loop

            if loop['loop_id'][0:2] == "J3":
                # try to find a good match quickly
                Q1 = create_query(loop, ifedata, discrepancy = 0.5)
                queries[loop['loop_id']]['Q1'] = Q1

                # try again
                Q2 = create_query(loop, ifedata, discrepancy = 0.8)
                queries[loop['loop_id']]['Q2'] = Q2
            elif loop['loop_id'][0] == "J":
                # first query puts > constraint in every position because homologous junctions
                Q1 = create_query(loop, ifedata, discrepancy = 3.0, increasing = True)
                queries[loop['loop_id']]['Q1'] = Q1

                # second query hopes to find a low-discrepancy match quickly
                Q2 = create_query(loop, ifedata, discrepancy = 0.5)
                queries[loop['loop_id']]['Q2'] = Q2

                if loop['loop_id'][1] in ["6","7","8","9"]:
                    # unrestricted order search with discrepancy 1 is too slow for these for now
                    Q = create_query(loop, ifedata, discrepancy = 1.0, increasing = True)
                    queries[loop['loop_id']]['Q'] = Q

    return queries, search_spaces, flanking_bp_queries, loops


def analyze_single_bulged_base_loops(candidates, query, search_space):
    not_rejected = []
    for candidate in candidates:
        DQ_code = []
        # query_bulge = get_bulge(query)
        # search_space_bulge = get_bulge(search_space)

        # if query_bulge != search_space_bulge:
        #     DQ_code.append(MISMATCHED_BULGE)

        filtered_candidate = {}

        if len(DQ_code) == 0:
            filtered_candidate['dq'] = []
            filtered_candidate['query_id'] = query['loop_info']['loop_id']
            filtered_candidate['search_space_id'] = search_space['loop_info']['loop_id']
            filtered_candidate['indices'] = candidate['indices']
            filtered_candidate['query_unit_ids'] = query['Q']['unitID']
            filtered_candidate['target_unit_ids'] = candidate['unitids']
            filtered_candidate['discrepancy'] = candidate['discrepancy']
            filtered_candidate['numStrands'] = query['Q']['numStrands']

            # filtered_candidate['query_rfam'] = query['loop_info']['rfam']
            # filtered_candidate['search_space_rfam'] = search_space['loop_info']['rfam']

        else: # if has a real DQ_code
            filtered_candidate['dq'] = list(set(DQ_code))
            filtered_candidate['discrepancy'] = candidate['discrepancy']
            filtered_candidate['numStrands'] = query['Q']['numStrands']

        not_rejected.append(filtered_candidate)

    temp_result = keep_lowest_discrepancy_candidate(not_rejected)
    return temp_result


def all_against_all_searches(loop_type,queries,search_spaces,flanking_bp_queries,load_saved_searches = True,reversed_search=False, save_path = DATAPATHRESULTS):
    """
    Search all loops against all loops
    Take in query strands and search space strands
    Flow:
    Set query_pbd_id to empty to check the first iteration
    Sort all pbd ids
    Load previous search results if query_pbd_id changes
    Save the results if search_pbd_id changes and new_results = true
    results is a dictionary store temporary search results
    search_results is the accumulated search results, get returned
    """

    new_results = False
    timer_data = myTimer("All against all search")
    search_space_pdb_id = ""
    search_results = {} # search_results store the final, accumulated results

    loop_counter = 0

    for search_space_id in sorted(search_spaces.keys(), reverse = reversed_search):
        loop_counter += 1

        search_space_length = len(search_spaces[search_space_id]['ifedata']['index_to_id'])
        search_space_bulged_length = len(search_spaces[search_space_id]['loop_info']['bulged'])
        is_single_bulged_search_space = search_space_bulged_length==1 and search_space_length==5

        start_time_on_this_pdb = time()

        # Only load previous result when there is a change in search_space_pdb_id.
        # Save the results if not the first iteration.
        if search_space_pdb_id != search_space_id.split("_")[1]:
            # save only once per pdb because that is faster
            if search_space_pdb_id != "":
                if new_results:
                    save_search_results_one_pdb(loop_type, search_space_pdb_id, results, save_path)
            search_space_pdb_id = search_space_id.split("_")[1]
            results = {} # results: {(query_id,search_space_id) : [dq,discrepancy]}
            if load_saved_searches:
                timer_data = myTimer("Loading files")
                # print('Loading search results for %s from %s' % (search_space_pdb_id, save_path))
                results = load_previous_search_results_one_pdb(loop_type, search_space_pdb_id, save_path)
            new_results = False

        print("Searching for %d query loops inside loop %s with size %2d, %4d/%4d, reporting hits only" % (len(queries),search_space_id,search_space_length,loop_counter,len(search_spaces)))
        for query_id in sorted(queries.keys()):

            if query_id != search_space_id:

                query_length = len(queries[query_id]['Q']['fullUnits'])
                query_bulge_length = len(queries[query_id]['loop_info']['bulged'])
                is_single_bulged_candidate = query_bulge_length==1 and query_length==5

                check_again = False

                # Temporary way to fix some missed flanking searches
                # When failed flanking it looks like [{'dq': [9], 'discrepancy': 99}]
                # if (query_id, search_space_id) in results:
                #     if type(results[(query_id, search_space_id)]) == list:
                #         if len(results[(query_id, search_space_id)]) > 0:
                #             if type(results[(query_id, search_space_id)][0]) == dict:
                #                 if "dq" in results[(query_id, search_space_id)][0]:
                #                     if results[(query_id, search_space_id)][0]["dq"] == [9]:
                #                         check_again = True
                #             elif results[(query_id, search_space_id)][0] == 9:
                #                 check_again = True
                #     elif results[(query_id, search_space_id)] == 9:
                #         check_again = True

                # re-check all single base bulge IL
                # if is_single_bulged_candidate:
                #     check_again = True

                # fix a problem where some loops had a difference between position and position_2023 in loop_positions
                # if queries[query_id]["loop_info"].get("check_again",False):
                #     check_again = True
                # if search_spaces[search_space_id]["loop_info"].get("check_again",False):
                #     check_again = True

                # revised HL criteria to insist on a good number of core positions
                # needs to be applied here to fix old queries that searched and found a match
                # if loop_type == 'HL' and query_length - query_bulge_length <= search_space_length:
                #     query_non_bulged_length = query_length - query_bulge_length
                #     search_space_non_bulged_length = search_space_length - search_space_bulged_length
                #     if (2 * query_non_bulged_length < search_space_non_bulged_length or 3 * query_non_bulged_length < search_space_length):
                #         print("query %s has %d non-bulge positions, search space %s has %d non-bulge, so skip" % (query_id,query_non_bulged_length,search_space_id,search_space_non_bulged_length))
                #         results[(query_id, search_space_id)] = SEARCH_SPACE_CONFLICT
                #         search_results[(query_id,search_space_id)] = [{'dq': [NO_CANDIDATES], 'discrepancy': 99}]
                #         continue

                # temporary: add more nucleotides to HL that have only 2 non-bulged positions, re-compute discrepancies
                # if loop_type == 'HL' and queries[query_id].get('added_nts',False):
                #     check_again = True

                # temporary: recalculate discrepancies of successful searches using flip
                # previous_discrepancy = 0
                # if (query_id, search_space_id) in results and (loop_type == 'J3' or (loop_type == 'IL' and search_space_length <= 20 and search_space_id < 'IL_8C3A_485')):
                #     if type(results[(query_id, search_space_id)][0]) == dict:
                #         check_again = True
                #         previous_discrepancy = results[(query_id, search_space_id)][0]['discrepancy']
                #         # print('Re-searching %s inside %s' % (query_id, search_space_id))

                # revised criterion so we don't lose core positions to an instance where a core position happens to be bulged
                # the instances may find each other when current search space becomes the query
                # needs to be applied here to fix old queries that searched and found a match
                # query_non_bulged_length = query_length - query_bulge_length
                # search_space_non_bulged_length = search_space_length - search_space_bulged_length
                # if query_non_bulged_length < search_space_non_bulged_length:
                #     # print("query %s has %d non-bulge positions, search space %s has %d non-bulge, so skip" % (query_id,query_non_bulged_length,search_space_id,search_space_non_bulged_length))
                #     results[(query_id, search_space_id)] = SEARCH_SPACE_CONFLICT
                #     search_results[(query_id,search_space_id)] = [{'dq': [NO_CANDIDATES], 'discrepancy': 99}]
                #     continue

                # for 3x3 IL, do not use discrepancy_flip, recompute those discrepancies
                # if queries[query_id]['Q'].get("3x3",False):
                #     check_again = True

                if (query_id, search_space_id) not in results or check_again:

                    # print('Searching for %s inside %s' % (query_id, search_space_id))
                    previous_discrepancy = 99
                    if queries[query_id]['Q'].get("3x3",False) and (query_id, search_space_id) in results:
                        previous_discrepancy = results[(query_id, search_space_id)][0]['discrepancy']

                    new_results = True

                    if loop_type == 'IL' and is_single_bulged_candidate:
                        # treat single base bulge IL differently, to keep them together according to bulged base

                        # default result
                        results[(query_id, search_space_id)] = MISMATCHED_BULGE
                        if is_single_bulged_search_space:

                            query_bulge = get_bulge(queries[query_id])
                            search_space_bulge = get_bulge(search_spaces[search_space_id])

                            if query_bulge == search_space_bulge:
                                # 5 nucleotide IL have a special query with discrepancy 98
                                # it uses all five nucleotides, so they have to be aligned
                                Q, candidates, elapsed_time = FR3D_search(Q = queries[query_id]['Q'],
                                            ifedata = search_spaces[search_space_id]['ifedata'], ifename = search_space_id,
                                            timerData = timer_data)
                                if candidates:
                                    temp_result = analyze_single_bulged_base_loops(candidates=candidates,query=queries[query_id],search_space=search_spaces[search_space_id])
                                    if temp_result:
                                        temp_result[0]['match_type'] = "geometric"
                                        results[(query_id, search_space_id)] = temp_result

                                    # print('SBL discrepancy %0.4f' % temp_result[0]['discrepancy'])
                                else:
                                    # I don't know why we would ever get here
                                    # print the query and search space loop ids
                                    print('Single bulge but no candidates found for %s inside %s' % (query_id, search_space_id))

                    elif search_space_id[0:2] =="IL" and is_single_bulged_search_space:
                        results[(query_id, search_space_id)] = MISMATCHED_BULGE

                    else:
                        query_non_bulged_length = query_length - query_bulge_length
                        search_space_non_bulged_length = search_space_length - search_space_bulged_length
                        if loop_type[0] == "J" and search_space_length - query_length > 15:
                            # pick the number above to cull out the slowest, most doomed searches
                            # but not much else
                            print('Skipping %s in %s because search space has %d nt and query has %d' % (query_id, search_space_id, search_space_length, query_length))
                            results[(query_id, search_space_id)] =  SEARCH_SPACE_CONFLICT
                        elif loop_type == "IL" and search_space_length > 35 and query_length < 20:
                            # avoid searching a not huge IL inside a huge IL
                            print('Skipping %s in %s because search space has %d nt and query has %d' % (query_id, search_space_id, search_space_length, query_length))
                            results[(query_id, search_space_id)] =  SEARCH_SPACE_CONFLICT
                        elif query_non_bulged_length < search_space_non_bulged_length:
                            # don't lose core positions to an instance where a core position happens to be bulged
                            results[(query_id, search_space_id)] = SEARCH_SPACE_CONFLICT
                        elif query_non_bulged_length > search_space_length:
                            # fixed positions of query are larger than search space
                            results[(query_id, search_space_id)] = SEARCH_SPACE_CONFLICT
                        elif  loop_type == 'HL' and (2 * query_non_bulged_length < search_space_non_bulged_length\
                            or 3 * query_non_bulged_length < search_space_length):
                            print('Skipping %s in %s because query has %d non-bulged and search space has %d non-bulged and size %d' % (query_id, search_space_id, query_non_bulged_length, search_space_non_bulged_length, search_space_length))
                            results[(query_id, search_space_id)] = SEARCH_SPACE_CONFLICT
                        else:
                            # Check flanking bp for IL and larger to make sure the match is plausible
                            if loop_type != 'HL':

                                # if query_id == 'J3_8VTW_036' and search_space_id == 'J3_5J7L_036':
                                #     print('About to do a flanking bp search J3_8VTW_036 inside J3_5J7L_036')
                                #     print(flanking_bp_queries[query_id]['Q'])

                                # if query_id == 'J3_5J7L_036' and search_space_id == 'J3_8VTW_036':
                                #     print('About to do a flanking bp search J3_5J7L_036 inside J3_8VTW_036')
                                #     print(queries[query_id]['Q'])

                                # print_dictionary(flanking_bp_queries[query_id]['Q']['interactionMatrix'])
                                # print("Flanking query above is for %s size %d inside %s size %d" % (query_id, query_non_bulged_length, search_space_id, search_space_length))

                                Q, candidates, elapsed_time = FR3D_search(Q=flanking_bp_queries[query_id]['Q'],ifedata = search_spaces[search_space_id]['ifedata'], ifename = search_space_id,
                                timerData = timer_data)

                                # if len(candidates) > 0:
                                #     print('Flanking found %s inside %s in %3d ways' % (query_id,search_space_id,len(candidates)))

                            if loop_type == "HL" or candidates:

                                # print('Searching for %s inside %s' % (query_id, search_space_id))

                                # do a full search for query inside of search_space
                                if loop_type[0] == "J":
                                    start_time = time()

                                    print_dictionary(queries[query_id]['Q1']['interactionMatrix'])
                                    print("Q1 query above is for %s size %d inside %s size %d" % (query_id, query_non_bulged_length, search_space_id, search_space_length))

                                    # do a quicker search with lower discrepancy first
                                    Q, candidates, elapsed_time = FR3D_search(Q = queries[query_id]['Q1'],
                                            ifedata = search_spaces[search_space_id]['ifedata'], ifename = search_space_id,
                                            timerData = timer_data)

                                    # print_dictionary(queries[query_id]['Q1'])
                                    print('Q1 found %d candidates in %9.4f seconds; query size %d, search space size %d' % (len(candidates),time() - start_time, Q['numPositions'], search_space_length))

                                    if len(candidates) == 0 and time() - start_time < 30:
                                        start_time = time()

                                        print_dictionary(queries[query_id]['Q2']['interactionMatrix'])
                                        print("Q2 query above is for %s size %d inside %s size %d" % (query_id, query_non_bulged_length, search_space_id, search_space_length))

                                        Q, candidates, elapsed_time = FR3D_search(Q = queries[query_id]['Q2'],
                                                ifedata = search_spaces[search_space_id]['ifedata'], ifename = search_space_id,
                                                timerData = timer_data)

                                        print('Q2 found %d candidates in %9.4f seconds; query size %d, search space size %d' % (len(candidates),time() - start_time, Q['numPositions'], search_space_length))

                                    if loop_type in ['J3','J4','J5'] and len(candidates) == 0 and time() - start_time < 400:
                                        start_time = time()

                                        print_dictionary(queries[query_id]['Q']['interactionMatrix'])
                                        print("Q query above is for %s size %d inside %s size %d" % (query_id, query_non_bulged_length, search_space_id, search_space_length))

                                        Q, candidates, elapsed_time = FR3D_search(Q = queries[query_id]['Q'],
                                                ifedata = search_spaces[search_space_id]['ifedata'], ifename = search_space_id,
                                                timerData = timer_data)

                                        print('Q  found %d candidates in %9.4f seconds; query size %d, search space size %d' % (len(candidates),time() - start_time, Q['numPositions'], search_space_length))

                                else:
                                    Q, candidates, elapsed_time = FR3D_search(Q = queries[query_id]['Q'],
                                            ifedata = search_spaces[search_space_id]['ifedata'], ifename = search_space_id,
                                            timerData = timer_data)

                                if candidates:
                                    temp_result = filter_out_conflicting_basepairs_and_stacks(candidates, queries[query_id], search_spaces[search_space_id])

                                    #print("%d candidates remain after conflicting basepairs" % len(temp_result))

                                    if temp_result:
                                        temp_result = filter_on_unmatched_nucleotides(temp_result,search_spaces[search_space_id],queries[query_id])

                                        # if len(temp_result) > 0:
                                        #     print("%d candidates remain after unmatched nucleotides" % len(temp_result))

                                        if temp_result:
                                            temp_result = keep_lowest_discrepancy_candidate(temp_result)
                                            dq = temp_result[0]['dq']
                                            if len(dq) == 0:
                                                temp_result[0]['match_type'] = "geometric"
                                                results[(query_id,search_space_id)] = temp_result

                                                if queries[query_id]['Q'].get("3x3",False):
                                                    new_discrepancy = temp_result[0]['discrepancy']
                                                    print("3x3 IL discrepancy for %s in %s goes from %0.2f to %0.2f" % (query_id,search_space_id,previous_discrepancy,new_discrepancy))

                                                # if check_again:
                                                #     new_discrepancy = temp_result[0]['discrepancy']
                                                #     if previous_discrepancy - new_discrepancy > 0.1:
                                                #         print("%s in %s discrepancy drops from %0.2f to %0.2f" % (query_id,search_space_id,previous_discrepancy,new_discrepancy))
                                            else:
                                                results[(query_id,search_space_id)] = dq[0]  # save one disqualification code
                                            # else:
                                            #     results[(query_id, search_space_id)] = UNMATCHED_BASEPAIR
                                        else:
                                            results[(query_id, search_space_id)] = UNMATCHED_BASEPAIR
                                    else:
                                        results[(query_id, search_space_id)] = CONFLICTING_BASEPAIRS_AND_STACKS

                                    if len(temp_result) > 0:
                                        print('Found %s inside %s in %3d ways, and %2d remain after filtering, %5d total results for %s' % (query_id,search_space_id,len(candidates),len(temp_result), len(results), search_space_pdb_id))

                                else: # no candidates from FR3D search
                                    #results[(query_id, search_space)] = [{'dq': 0, 'discrepancy': 99}] # hoping disqualification code of 0 makes sense for NO MATCH
                                    #Put 0 instead of list
                                    results[(query_id, search_space_id)] = NO_CANDIDATES
                                # after candidate sorting, we will no longer use Find_lowest()
                                # run test_motif_atlas_code::check_interaction() after EACH FR3D search
                            else:
                                results[(query_id, search_space_id)] =  FLANKING_BP_CONFLICT

                if isinstance(results[(query_id, search_space_id)],int):
                    search_results[(query_id,search_space_id)] = [{'dq': [NO_CANDIDATES], 'discrepancy': 99}]
                else:
                    search_results[(query_id,search_space_id)] = results[(query_id,search_space_id)]


            if time() - start_time_on_this_pdb > 120:
                save_search_results_one_pdb(loop_type, search_space_pdb_id, results, save_path)
                print('Saved %d results for %s and %s in %s' % (len(results), loop_type, search_space_pdb_id, save_path))
                start_time_on_this_pdb = time()


    # save the results for the last PDB, which may be the only PDB
    if len(results) > 0:
        print('Saving %d results for %s and %s in %s' % (len(results), loop_type, search_space_pdb_id, save_path))
        save_search_results_one_pdb(loop_type, search_space_pdb_id, results, save_path)

    return search_results


def clique_criterion(clique,disc_m,cluster_method):
    """
    Apply a criterion to evaluate a clique
    """

    if 'average' in cluster_method:
        sum_disc = 0
        counter = 0.0001   # force decimal division, avoid division by zero
        for i in clique:
            for j in clique:
                if j > i:
                    sum_disc += disc_m[i][j]
                    counter += 1

        return sum_disc/counter

    elif 'complete' in cluster_method:
        max_disc = 0
        for i in clique:
            for j in clique:
                if j > i:
                    max_disc = max(max_disc,disc_m[i][j])
        return max_disc

    return 0


def find_best_clique(disc_m,clique_list,ratio,cluster_method,descending=True):
    if not descending:
        clique_list.sort(key=len,reverse=True)

    # set a size limit determined by the ratio
    size_limit = round(len(clique_list[0])*ratio)

    # largest clique has the first opportunity to be chosen
    new_clique = clique_list[0]
    min_criterion = clique_criterion(new_clique,disc_m,cluster_method)

    counter = 0
    for clique in clique_list:
        if len(clique) >= size_limit and len(clique) > 1:
            counter += 1
            disc = clique_criterion(clique,disc_m,cluster_method)
            if disc < min_criterion:
                #print(f"Switching to clique of size {len(i)} with disc {disc} over {min_criterion}")
                min_criterion = disc
                new_clique = clique
        else:
            break
    print("Found clique of size %3d with criterion %4.2f going over %3d cliques from size %3d down to %3d out of %4d cliques" %
        (len(new_clique), min_criterion, counter, len(clique_list[0]), size_limit, len(clique_list)))

    return new_clique


def cluster_motifs3(clique_list,disc_m,ratio,cluster_method):
    """
    Select the best large clique, remove it, select the next best, continue.
    clique_list is a list of all largest cliques, each of which is a set
    """

    clique_list = [set(clique) for clique in clique_list]

    groups = []
    counter = 0

    while True:
        # sort from largest clique to smallest
        clique_list.sort(key=len,reverse=True)
        if len(clique_list[0]) == 1 or len(clique_list) == 0 or len(clique_list[0]) == 0:
            break
        best_clique = find_best_clique(disc_m,clique_list,ratio,cluster_method)
        groups.append(list(best_clique))
        counter += len(best_clique)

        # remove loops in best_clique from all other cliques
        #clique_list = [[x for x in clique if x not in best_clique] for clique in clique_list]
        clique_list = [clique - best_clique for clique in clique_list]
        # remove empty cliques, if any
        clique_list = [clique for clique in clique_list if len(clique) > 0]

    # add in singleton cliques that remain, each one only once
    singletons = set()
    for clique in clique_list:
        if len(clique) == 1:
            singletons.add(list(clique)[0])

    print('%d non-singleton groups have a total of %4d loops' % (len(groups),counter))
    print("Found %3d singletons" % len(singletons))

    for singleton in singletons:
        groups.append([singleton])

    return groups


def cluster_loops_maximum_clique(MM,disc_m,dq_matrix,ratio=1,cluster_method='clique_average'):
    """
    Code from Seyoung.
    MM is a matrix of 0 and 1, 0 means discrepancy > 1, 1 means discrepancy < 1.
    Find the largest clique of size S
    Compare to cliques down to size ratio*S to find the tightest clique, then remove it
    Repeat until no cliques are left
    """

    # get a list of all largest possible cliques in the matching matrix MM
    G = nx.from_numpy_matrix(MM)
    clique_list = list(nx.find_cliques(G))

    # choose the tightest cliques from largest to smallest
    groups = cluster_motifs3(clique_list,disc_m,ratio,cluster_method)

    #heat_map(cliques,groups,disc_m)
    # if show_dq_codes:
    #     get_dq_between_motif_groups(cliques,dq_matrix,groups)

    return groups


def calculate_distance_matrices(search_results, all_loop_ids):
    """
    Calculate matching matrix, discrepancy matrix, disqualification matrix
    """

    # construct a "matching matrix" which has values 0 and 1, 1 for matches
    N = len(all_loop_ids)

    MM = np.zeros(shape=(N,N))       # matching matrix
    disc_m = np.zeros(shape=(N,N))   # discrepancy matrix
    dq_matrix = [[-1] * N for i in range(N)]

    for i in range(0,N):
        loop_i_id = all_loop_ids[i]
        for j in range(i+1,N):
            loop_j_id = all_loop_ids[j]

            i_j_disqualification = search_results[(loop_i_id,loop_j_id)][0]['dq']
            j_i_disqualification = search_results[(loop_j_id,loop_i_id)][0]['dq']

            if len(i_j_disqualification) == 0 or len(j_i_disqualification) == 0:
                # the loops can be in the same motif group
                MM[i][j] = 1
                i_j_discrepancy = search_results[(loop_i_id,loop_j_id)][0]['discrepancy']
                j_i_discrepancy = search_results[(loop_j_id,loop_i_id)][0]['discrepancy']
                disc_m[i][j] = min(i_j_discrepancy,j_i_discrepancy)
            else:
                # the loops cannot be in the same motif group
                MM[i][j] = 0
                dq_matrix[i][j] = i_j_disqualification[0]
                dq_matrix[j][i] = j_i_disqualification[0]
                disc_m[i][j] = 99999999999.0                # force them apart

            MM[j][i] = MM[i][j]
            disc_m[j][i] = disc_m[i][j]

    return disc_m, MM, dq_matrix


def cluster_loops_hierarchical(disc_m,cluster_method,all_loop_ids):
    """
    Use hierarchical clustering to find tight groups of loops
    Output is a list of motif groups, and each motif group is a
    list of indices into all_loop_ids.
    """

    N = len(disc_m)

    method_name = cluster_method.replace("hierarchical_","")

    # get a list of mergers between individuals and clusters
    Z = linkage(squareform(disc_m), method_name)

    # track the different groups as they are merged together
    group_num_to_list = {}   # the members of the group, but emptied once merged

    for i in range(0,N):
        # the first groups are singletons
        group_num_to_list[i] = [i]

    # loop over the mergers that would happen if you followed hierarchical clustering
    # all the way to the end, and construct the groups that would be formed
    for merger in Z:
        i += 1
        if merger[2] < 10:
            # single bulged bases can have discrepancy above 1
            # just avoid whatever large number indicates a disqualification
            # once merger[2] gets above 10, any further mergers would be higher
            # than 10 as well, so we never see them

            a = int(merger[0])
            b = int(merger[1])

            # if the merger that was to construct a was skipped, don't merge further
            if len(group_num_to_list[a]) == 0:
                group_num_to_list[i] = []
                continue
            if len(group_num_to_list[b]) == 0:
                group_num_to_list[i] = []
                continue

            min_disc_a_b = 99999999999.0
            for ii in group_num_to_list[a]:
                for jj in group_num_to_list[b]:
                    min_disc_a_b = min(min_disc_a_b,disc_m[ii][jj])

            print('Considering a merger at %s cluster distance %0.4f and min distance %0.4f between group of size %d and of size %d' % (method_name,merger[2],min_disc_a_b,len(group_num_to_list[a]),len(group_num_to_list[b])))

            max_disc = 0
            for ii in group_num_to_list[a]:
                for jj in group_num_to_list[a]:
                    max_disc = max(max_disc,disc_m[ii][jj])
            max_disc_a_a = max_disc

            max_disc = 0
            for ii in group_num_to_list[b]:
                for jj in group_num_to_list[b]:
                    max_disc = max(max_disc,disc_m[ii][jj])
            max_disc_b_b = max_disc

            loop_ids = ['A max_disc_a_a %0.4f' % max_disc_a_a]
            for ww in group_num_to_list[a]:
                loop_ids.append(all_loop_ids[ww])
            print(" ".join(sorted(loop_ids)))
            loop_ids = ['A max_disc_b_b %0.4f' % max_disc_b_b]
            for ww in group_num_to_list[b]:
                loop_ids.append(all_loop_ids[ww])
            print(" ".join(sorted(loop_ids)))

            if len(group_num_to_list[a]) >= 3 and len(group_num_to_list[b]) >= 3:
                if merger[2] > 0.6 and min_disc_a_b > 0.6 and abs(merger[2]-min_disc_a_b) < 0.2:
                    if max_disc_a_a < 0.4 and max_disc_b_b < 0.4:
                        # groups that deserve to be split and satisfy these criteria
                        # https://rna.bgsu.edu/rna3dhub/motif/view/IL_93341.4
                        # Considering a merger at average cluster distance 0.7165 and min distance 0.6817 between group of size 4 and of size 3
                        # A max_disc_a_a 0.1341 IL_5T83_002 IL_5U3G_002 IL_6DLR_002 IL_7MLW_003
                        # A max_disc_b_b 0.2188 IL_5T83_006 IL_6DLR_006 IL_7MLW_009
                        # https://rna.bgsu.edu/rna3dhub/motif/view/IL_48283.1
                        # Considering a merger at average cluster distance 0.7060 and min distance 0.6656 between group of size 6 and of size 8
                        # A max_disc_a_a 0.1607 IL_4V9F_063 IL_4WF9_066 IL_5J7L_308 IL_7A0S_063 IL_8P9A_307 IL_9DFE_068
                        # A max_disc_b_b 0.2539 IL_1NBS_007 IL_4V9F_007 IL_4V9F_049 IL_4V9F_092 IL_4V9F_106 IL_5J7L_251 IL_8P9A_289 IL_9DFE_007
                        # https://rna.bgsu.edu/rna3dhub/motif/view/IL_31461.3
                        # Considering a merger at average cluster distance 0.7062 and min distance 0.6415 between group of size 3 and of size 5
                        # A max_disc_a_a 0.1910 IL_4LFB_013 IL_5J7L_016 IL_6CZR_125
                        # A max_disc_b_b 0.1858 IL_4LFB_015 IL_5J7L_018 IL_6CZR_127 IL_8C3A_416 IL_8P9A_397

                        # groups that might not deserve to be split
                        # splitting these two prevents both of them from being merged later into
                        # a group they belong in
                        # https://rna.bgsu.edu/rna3dhub/motif/view/IL_31462.3
                        # Considering a merger at average cluster distance 0.8985 and min distance 0.7292 between group of size 5 and of size 7
                        # A max_disc_a_a 0.4259 IL_4V9F_101 IL_4WF9_112 IL_5J7L_355 IL_7DLZ_002 IL_9DFE_110
                        # A max_disc_b_b 0.8059 IL_4LFB_002 IL_4P97_002 IL_4PHY_002 IL_5CNR_002 IL_5J7L_002 IL_6CZR_114 IL_9DFE_100

                        # groups that might not deserve to be split

                        print('Skipping this merger by making the merged group empty')
                        group_num_to_list[i] = []
                        continue

            group_num_to_list[i] = group_num_to_list[a] + group_num_to_list[b]
            # now that these are merged, empty them so we don't put one loop into
            # multiple groups
            group_num_to_list[a] = []
            group_num_to_list[b] = []

    groups = [g for g in group_num_to_list.values() if len(g) > 0]

    groups = sorted(groups, key = lambda x : len(x), reverse = True)

    for i,g in enumerate(groups):
        if len(g) > 1:
            print('Group %3d has size %3d and includes %s' % (i,len(g),g[0]))

    print('Found %4d groups' % len(groups))

    return groups


def order_group_by_similarity(cluster,disc_m):
    """
    Re-order each of the motif clusters by similarity
    """

    if len(cluster) > 1:
        dist_matrix = np.zeros(shape=(len(cluster),len(cluster)))

        # get the relevant part of the distance matrix
        for i in range(0,len(cluster)):
            loop_i = cluster[i]
            for j in range(i,len(cluster)):
                loop_j=cluster[j]

                dist_matrix[i][j] = disc_m[loop_i][loop_j]
                dist_matrix[j][i] = disc_m[loop_i][loop_j]

        order = treePenalizedPathLength(dist_matrix)

        reordered_dist_matrix = np.zeros(shape=(len(cluster),len(cluster)))
        for i in range(0,len(cluster)):
            for j in range(0,len(cluster)):
                reordered_dist_matrix[i][j] = dist_matrix[order[i]][order[j]]

        reordered_cluster = [cluster[x] for x in order]
    else:
        reordered_cluster = cluster
        reordered_dist_matrix = np.zeros(shape=(1,1))

    return reordered_cluster, reordered_dist_matrix


def search_for_centroid_loop(loop_ids,loop_ids_to_discrepancy,search_results,queries):
    """
    Identify the centroid loop in a motif group
    Return the loop_id of the centroid
    Also return
    """

    min_total_disc = maxsize
    centroid_id = ""
    centroid_units = []
    for oneLoop in loop_ids:
        oneLoopUnits = set(queries[oneLoop]["Q"]["unitID"]) # unit ids used for query, not bulged
        total_disc = 0
        largestNumberOfAlignment = 0
        for anotherLoop in loop_ids:
            total_disc += loop_ids_to_discrepancy[oneLoop][anotherLoop]
            if oneLoop != anotherLoop:
                # oneLoop is query_id, and anotherLoop is search_space_id
                # query_unit_ids belong with query_id, target_unit_ids belong with search_space_id
                # oneLoop -- query_id --query_unit_ids
                # anotherLoop --- search_space_id --- target_unit_ids
                sr = search_results[(oneLoop,anotherLoop)][0]

                if "target_unit_ids" in sr:
                    # print('This is the first sr, %s \n the oneLoop is %s' % (sr, oneLoop))
                    srOneLoopUnits = sr["query_unit_ids"] # check if its from oneLoop: query_unit_ids belongs oneLoop, target_unit_ids belongs anotherLoop
                    oneLoopUnits = oneLoopUnits.intersection(set(srOneLoopUnits)) # intersection between query unit ids and oneLoopUnits

                else:
                    sr = search_results[(anotherLoop,oneLoop)][0]
                    # print('This is the second sr, %s \n the oneLoop is %s' % (sr, oneLoop))
                    srOneLoopUnits = sr["target_unit_ids"] # oneLoop
                    oneLoopUnits = oneLoopUnits.intersection(set(srOneLoopUnits)) # same thing down here

        # keep the loop_id of the loop with more aligned nts and smaller total/average discrepancy to all others
        if len(oneLoopUnits) >= largestNumberOfAlignment and (total_disc < min_total_disc):
            min_total_disc = total_disc
            largestNumberOfAlignment = len(oneLoopUnits)
            centroid_id = oneLoop
            # store the units in the original order, not the random order after set intersection
            centroid_units = []
            for unit_id in queries[oneLoop]["Q"]["unitID"]:
                if unit_id in oneLoopUnits:
                    centroid_units.append(unit_id)

    return centroid_id, centroid_units


# def align_candidate_units_to_centroid_units(centroid_core_unit_ids, centroid_unit_ids, candidate_unit_ids, centroid_to_candidate, core_positions_candidate_id):
#     for centroid_unit_id in centroid_core_unit_ids:
#         centroid_unit_id_index = centroid_unit_ids.index(centroid_unit_id)
#         candidate_unit_id = candidate_unit_ids[centroid_unit_id_index]
#         centroid_to_candidate[centroid_unit_id] = candidate_unit_id
#         core_positions_candidate_id.append(centroid_to_candidate[centroid_unit_id])


def align_nts(loop_ids, loop_ids_to_discrepancy, search_results, queries):
    """
    Identify centroid instance in the given cluster and its core nucleotides and align all instances to that
    loop_ids is a list of text strings
    loop_ids_to_discrepancy[loop_1_id][loop_2_id] = discrepancy
    """

    # identify centroid loop id, and get its units in the right order
    centroid_id, centroid_units = search_for_centroid_loop(loop_ids, loop_ids_to_discrepancy, search_results, queries)

    print('check order:',centroid_id, centroid_units)

    # start to count how long each strand is, in the centroid
    dict_for_separate_strand = {}
    for strand_index in range(len(queries[centroid_id]["loop_info"]["strand"])):
        dict_for_separate_strand[strand_index] = []

    # need OneLoopUnits around here

    # add a condition for skipping HLs ? really need it ?
    # for strand in queries[centroid_id]["loop_info"]["strand"]:
    #     strand_name = queries[centroid_id]["loop_info"]["strand"].index(strand)
    #     for unit_id in strand:
    #         if unit_id in centroid_units:
    #             dict_for_separate_strand[strand_name].append(unit_id)

    # add non-bulged unit ids to each strand; dict_for_unit_sets omits bulged nucleotides
    for unit_id in centroid_units:
        for strand_index, strand in enumerate(queries[centroid_id]["loop_info"]["strand"]):
            if unit_id in strand:
                dict_for_separate_strand[strand_index].append(unit_id)

    # find which strand is the longest, accounting for bulged nucleotides
    max_key = max(dict_for_separate_strand, key=lambda k: len(dict_for_separate_strand[k]))
    centroid_core_unit_ids = []
    if max_key != 0 and centroid_id.startswith('IL'):
        # longest strand is not the first strand
        print('Changing order of IL strands to put the longest one first')
        original_order = list(dict_for_separate_strand.keys())
        new_order = original_order[max_key:] + original_order[:max_key]
        sorted_dict_for_separate_strand = OrderedDict((k, dict_for_separate_strand[k]) for k in new_order)
        for strand_name, unit_ids in sorted_dict_for_separate_strand.items():
            centroid_core_unit_ids.extend(unit_ids)
    else:
        # original_order = list(dict_for_separate_strand.keys())
        # new_order = original_order[max_key:] + original_order[:max_key]
        for strand_name, unit_ids in dict_for_separate_strand.items():
            centroid_core_unit_ids.extend(unit_ids)

    loop_id_to_core_units = {}
    loop_id_to_core_units[centroid_id] = centroid_core_unit_ids

    # align other loops with centroid loop
    for candidate_id in loop_ids:
        if candidate_id == centroid_id:
            continue

        sr = search_results[(centroid_id,candidate_id)][0]
        if "target_unit_ids" in sr:
            candidate_unit_ids = sr["target_unit_ids"]
            centroid_unit_ids = sr["query_unit_ids"]
        else:
            sr = search_results[(candidate_id,centroid_id)][0]
            candidate_unit_ids = sr["query_unit_ids"]
            centroid_unit_ids = sr["target_unit_ids"]

        # align_candidate_units_to_centroid_units(centroid_core_unit_ids, centroid_unit_ids, candidate_unit_ids, centroid_to_candidate, loop_id_to_core_units[candidate_id])

        loop_id_to_core_units[candidate_id] = []
        for centroid_unit_id in centroid_core_unit_ids:
            # because strands can change order, index in centroid_core_unit_ids may differ from in centroid_unit_id
            centroid_unit_id_index = centroid_unit_ids.index(centroid_unit_id)
            candidate_unit_id = candidate_unit_ids[centroid_unit_id_index]
            loop_id_to_core_units[candidate_id].append(candidate_unit_id)

    return loop_id_to_core_units, centroid_id


def evaluate_grouping(motif_groups,loop_id_to_annotation):
    """
    Give basic statistics about the group sizes
    For each loop annotation, tell the largest number of times it occurs in one group
    For each loop annotation, tell the number of groups it appears in
    """

    annotation_to_group_to_count = {}
    count_singletons = 0

    for i, motif_group in enumerate(motif_groups):
        if len(motif_group) == 1:
            count_singletons += 1

        for loop_id in motif_group:
            if loop_id in loop_id_to_annotation:
                annotation = loop_id_to_annotation[loop_id]

                if not annotation in annotation_to_group_to_count:
                    annotation_to_group_to_count[annotation] = {}

                if not i in annotation_to_group_to_count[annotation]:
                    annotation_to_group_to_count[annotation][i] = 0

                annotation_to_group_to_count[annotation][i] += 1


    output = []

    output.append("%d groups" % (len(motif_groups)))
    output.append("%d singletons" % count_singletons)
    output.append("TotalCount NumberOfGroups MaxCount Annotation")

    for annotation in sorted(annotation_to_group_to_count.keys()):
        total_count = sum(annotation_to_group_to_count[annotation].values())
        num_groups = len(annotation_to_group_to_count[annotation].keys())
        max_count = max(annotation_to_group_to_count[annotation].values())
        output.append("%d\t%d\t%d\t%s" % (total_count,num_groups,max_count,annotation))

    print("\n".join(output))

def keep_groups_with_x_ray_instance(motif_groups, all_loop_ids, pdb_list):
    """
    Only keep motif groups that have at least one instance from an x-ray structure

    motif_groups is a list of lists, each list is a list of indices of loop ids
    """

    keep_motif_groups = []
    for motif_group in motif_groups:
        keep_group = False
        for loop_index in motif_group:
            loop_id = all_loop_ids[loop_index]
            pdb_id = loop_id.split("_")[1]  # get the PDB id from the loop id
            if pdb_id in pdb_list:
                keep_group = True
                break
        if keep_group:
            keep_motif_groups.append(motif_group)
        else:
            print("Not keeping motif group with these loop ids: %s" % [all_loop_ids[x] for x in motif_group])

    return keep_motif_groups


def cluster_loops(loop_position_to_border_unit_id, output_dir = './', molecule_type = 'RNA', x_ray_pdb_list = []):
    """
    output_dir is like /usr/local/pipeline/hub-core/MotifAtlas/Releases/HL_3.87_2024-08-02_20:30
    """

    queries = {}
    search_spaces ={}
    flanking_bp_queries = {}

    # to search faster, run two instances of this program with True/False on next line
    reversed_search = False
    reversed_search = True

    cluster_method = 'clique_maximum'
    cluster_method = 'clique_average'
    cluster_method = 'hierarchical_complete'
    cluster_method = 'hierarchical_average'

    load_saved_searches = True   # load saved search results for each PDB file; faster!

    ratio = 0.9          # the ratio that determines what range of clique sizes to compare

    use_loop_annotations = True  # a diagnostic technique

    timer_data = myTimer("start")
    timer_data = myTimer('Set up loops')

    loop_ids = sorted(loop_position_to_border_unit_id.keys())
    loop_type = loop_ids[0].split("_")[0]

    print("Grouping %d %s loops of interest" % (len(loop_ids),loop_type))

    # turn the raw list of loops and their nucleotides into a list of loop dictionaries
    loops = identify_strands(loop_position_to_border_unit_id)

    # sort the list of loops by loop id
    loops = sorted(loops, key = lambda x : x['loop_id'])

    timer_data = myTimer("Set up queries and ifedata",timer_data)

    # for each loop, set it up as a query Q and as a search space ifedata
    # also add bulged nucleotides as long as we are already loading pairwise interactions
    print('Creating queries and search spaces for %d loops' % len(loops))
    queries, search_spaces, flanking_bp_queries, loops = create_queries_and_search_spaces(loops)

    # index the loop ids in a list, so you can refer to them by number
    all_loop_ids = list(queries.keys())

    print("Found %d loops with full coordinate data" % len(all_loop_ids))

    print("Start all against all searches")
    timer_data = myTimer("All against all searches",timer_data)
    motif_result_path = DATAPATHATLAS
    search_results = all_against_all_searches(loop_type,queries,search_spaces,flanking_bp_queries,load_saved_searches,reversed_search, motif_result_path)

    timer_data = myTimer("Calculate distance matrices",timer_data)
    print("Calculating distance matrices")

    disc_m, MM, dq_matrix = calculate_distance_matrices(search_results, all_loop_ids)

    # diagnostic: load text annotations of loops, for later diagnostics
    if use_loop_annotations:
        # load loop annotations for testing
        loop_id_to_annotation = defaultdict(str)
        annotation_filename = os.path.join(output_dir, "loop_annotations.txt")
        if os.path.exists(annotation_filename):
            with open(annotation_filename,"r") as file:
                for line in file:
                    data = line.split("\t")
                    # map loop id to the first annotation; skip the broader annotation and author
                    loop_id_to_annotation[data[0]] = data[1]

    # diagnostic: put literally all of the loops into one ordering; can be a useful diagnostic
    order_all_loops = False
    if order_all_loops:
        timer_data = myTimer("Order all loops",timer_data)
        print("Ordering all loops")

        # sort loop instances by similarity and list them to the screen
        order = treePenalizedPathLength(disc_m,100)  # map new position to loop number
        for i in range(0,len(order)):
            j = order[i]
            loop_id = all_loop_ids[j]

            if i > 0:
                k = order[i-1]  # previous loop number
                d = disc_m[j][k]
            else:
                d = 0

            if d > 100:
                d = 99

            mg = loop_id_to_motif_id[loop_id]
            ann = loop_id_to_annotation[loop_id]

            print("%d\t%s\t%8.4f\t%s\t%s" % (i,loop_id,d,mg,ann))

    # diagnostic: compare to Matlab motif groups
    if False:
        # map loop id to index
        loop_id_to_index = {}
        for i, loop_id in enumerate(all_loop_ids):
            loop_id_to_index[loop_id] = i

        not_found = set()
        for motif_group, loop_ids in motif_id_to_loops.items():
            print('Motif group %s' % motif_group)
            for loop_id1 in loop_ids:
                if not loop_id1 in loop_id_to_index:
                    not_found.add(loop_id1)
                    continue
                i = loop_id_to_index[loop_id1]
                for loop_id2 in loop_ids:
                    if not loop_id2 in loop_id_to_index:
                        not_found.add(loop_id2)
                        continue
                    j = loop_id_to_index[loop_id2]

                    """
                    if disc_m[i][j] > 2:
                        print("  %s and %s discrepancy %8.4f" % (loop_id1,loop_id2,disc_m[i][j]))
                        print(search_results[(loop_id1,loop_id2)])
                    """

        # list loops whose coordinates are not available, if any
        for i, loop_id in enumerate(list(not_found)):
            print('  %s is not found here, #%d of %d' % (loop_id,i,len(not_found)))

    timer_data = myTimer("Cluster loops",timer_data)
    print("Clustering loops using %s" % cluster_method)

    if cluster_method in ['clique_average','clique_maximum']:
        motif_groups = cluster_loops_maximum_clique(MM,disc_m,dq_matrix,ratio,cluster_method)
    elif cluster_method in ['hierarchical_average','hierarchical_complete']:
        motif_groups = cluster_loops_hierarchical(disc_m,cluster_method,all_loop_ids)

    # optionally remove motif groups made only of cryo-EM structures
    if len(x_ray_pdb_list) > 0:
        motif_groups = keep_groups_with_x_ray_instance(motif_groups,all_loop_ids,x_ray_pdb_list)

    # diagnostic: find maximum distance within each group
    group_max_distance = {}
    for i, motif_group in enumerate(motif_groups):
        maxd = 0
        for a in motif_group:
            for b in motif_group:
                maxd = max(maxd,disc_m[a][b])
        group_max_distance[i] = maxd

    # diagnostic: find minimum distance between groups
    group_min_distance = {}
    for i, motif_group1 in enumerate(motif_groups):
        group_min_distance[i] = []
        for j, motif_group2 in enumerate(motif_groups):
            if i == j:
                continue
            mind = 1000
            for a in motif_group1:
                for b in motif_group2:
                    mind = min(mind,disc_m[a][b])
            if mind < 10:
                group_min_distance[i].append((mind,j,len(motif_group2)))
        group_min_distance[i] = sorted(group_min_distance[i])

        if len(motif_group1) > 1:
            print("Group %03d has %3d instances and intra-group max distance %6.2f and minimum distance" % (i,len(motif_group1),group_max_distance[i]))
            for mind,j,length in group_min_distance[i]:
                if mind < group_max_distance[i]:
                    print("  %6.2f to group %03d with size %3d" % (mind,j,length))

    # form motif groups by loop ids
    group_number_to_loop_id = {}
    loop_id_to_group_number = {}
    motif_groups_by_loop_id = []
    for i, motif_group in enumerate(motif_groups):
        group_number_to_loop_id[i] = [all_loop_ids[x] for x in motif_group]
        motif_groups_by_loop_id.append(group_number_to_loop_id[i])
        for x in motif_group:
            if all_loop_ids[x] in loop_id_to_group_number:
                print("%s was already assigned to group %3d, now it is being assigned to %3d" % (all_loop_ids[x],loop_id_to_group_number[all_loop_ids[x]],i))
            loop_id_to_group_number[all_loop_ids[x]] = i

    # diagnostic: evaluate how well the new clusters do, with loop annotations zzz
    if use_loop_annotations:
        print("Evaluation of how well the new clusters work")
        evaluate_grouping(motif_groups_by_loop_id,loop_id_to_annotation)

    # diagnostic: list loops that have low discrepancy but are in different groups, with loop annotations if any
    if use_loop_annotations:
        print("Loops with low discrepancy but in different groups")
        print("id1   #nts    group_index group_size  annotation1 discrepancy    id2 #nts    group2_index    group2_size")
        # HL_3IGI_008	6	2	56	UNCG	0.3946	HL_2Z75_004	6	0	305	GNRA

        similar_pairs = []
        for i, motif_group1 in enumerate(motif_groups):
            for j, motif_group2 in enumerate(motif_groups):
                if i <= j:
                    continue
                for a in motif_group1:
                    a_annotation = loop_id_to_annotation[all_loop_ids[a]]
                    a_length = len(loop_position_to_border_unit_id[all_loop_ids[a]])
                    for b in motif_group2:
                        b_annotation = loop_id_to_annotation[all_loop_ids[b]]
                        b_length = len(loop_position_to_border_unit_id[all_loop_ids[b]])
                        if disc_m[a][b] < 0.4:
                            t = "%s\t%d\t%d\t%d\t%20s\t%0.4f\t%s\t%d\t%d\t%d\t%s" % (all_loop_ids[a],a_length,i,len(motif_group1),a_annotation,disc_m[a][b],all_loop_ids[b],b_length,j,len(motif_group2),b_annotation)
                            similar_pairs.append(t)
                            print(t)
        with open(os.path.join(output_dir,loop_type+"_similar_pairs.txt"),"w") as file:
            file.write("\n".join(similar_pairs))

    # diagnostic: Did any loops not get mapped to a motif group?
    print("Mapped %4d loop ids to motif groups, out of %4d loop ids" % (len(loop_id_to_group_number.keys()),len(all_loop_ids)))
    missed_loop_ids = set(all_loop_ids) - set(loop_id_to_group_number.keys())
    if len(missed_loop_ids) > 0:
        print('There are %3d missing loop ids:' % len(missed_loop_ids))
        print(sorted(missed_loop_ids))

    # diagnostic: identify groups a loop could be in (including its own) and get max and mean discrepancy
    loop_id_to_matching_groups = defaultdict(dict)
    for i in range(0,len(all_loop_ids)):
        loop_id = all_loop_ids[i]
        # group_num_to_discrepancy = defaultdict(list)

        for group_num in range(0,len(motif_groups)):
            discrepancies = []
            for j in motif_groups[group_num]:
                dij = disc_m[i][j]
                if dij < 999:
                    discrepancies.append(dij)
                else:
                    discrepancies = []
                    break

            if len(discrepancies) > 0:
                loop_id_to_matching_groups[loop_id][group_num] = (max(discrepancies),np.mean(discrepancies))

    timer_data = myTimer("Write output",timer_data)
    print("Writing output")

    folder_name = output_dir + '/2ds'
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    # re-order the instances in each group, find core positions, write HTML file for inspection
    for num, motif_group in enumerate(motif_groups):

        group_num = num + 1

        # re-order distance matrices and change indexing from numeric to loop ids
        reordered_motif_group, reordered_dist_matrix = order_group_by_similarity(motif_group,disc_m)

        # index motif group and discrepancy matrix by loop ids
        reordered_loop_ids = [all_loop_ids[x] for x in reordered_motif_group]

        loop_ids_to_discrepancy = {}
        for i, loop_i_id in enumerate(reordered_loop_ids):
            loop_ids_to_discrepancy[loop_i_id] = {}
            for j, loop_j_id  in enumerate(reordered_loop_ids):
                loop_ids_to_discrepancy[loop_i_id][loop_j_id] = reordered_dist_matrix[i][j]

        # show_discrepancy(search_results, ['IL_3IGI_011'],['IL_4V9F_007'])

        # identify core positions in each group
        # Why?  The instances in a cluster might not all have the same number of nucleotides
        # It's common for a loop in one structure, one organism to have one or more "bulged out" nucleotides
        # The FR3D search technique we use is robust to those bulged out nucleotides
        # The FR3D search basically gives pairwise alignments between loop instances
        # But for a motif group, we need a complete, single alignment, not a bunch of pairwise alignments
        # We need to identify "core" and "non-core" nucleotides

        loop_id_to_core_units, centroid_id = align_nts(reordered_loop_ids, loop_ids_to_discrepancy, search_results, queries)

        print('Group %03d centroid %s core units %s' % (group_num, centroid_id, loop_id_to_core_units[centroid_id]))

        # collect interactions made in each core position
        # loop_id_type_pair_to_interaction[loop_id]['basepair'][corePositionPair] = 'cWW'
        # name could be loop_id_type_pair_to_interaction
        loop_id_type_pair_to_interaction = get_matrix_for_consensus_interactions(loop_id_to_core_units)

        # insist on cWW pairs in flanking positions in the loop
        centroid_unit_ids = loop_id_to_core_units[centroid_id]
        if loop_type == 'HL':
            # make a list using the first and last entry
            centroid_flanking_unit_ids = centroid_unit_ids[0:1] + centroid_unit_ids[-1:]
        else:
            centroid_flanking_unit_ids = flanking_bp_queries[centroid_id]["Q"]["unitID"]

        i = 0
        while i < len(centroid_unit_ids)-1:
            unit_id = centroid_unit_ids[i]
            if unit_id in centroid_flanking_unit_ids:
                if i == 0:
                    pair = "%d - %d" % (1,len(centroid_unit_ids))
                else:
                    pair = "%d - %d" % (i+1,i+2)
                    i += 1
                current_interaction = loop_id_type_pair_to_interaction[centroid_id]['basepair'].get(pair,'nul')
                print("%s flanking pair %s is %s should be cWW" % (centroid_id,pair,current_interaction))
                # save that pair for all instances, not just the centroid, then it gets through to the signature
                for loop_id in loop_id_type_pair_to_interaction.keys():
                    loop_id_type_pair_to_interaction[loop_id]['basepair'][pair] = 'cWW'
                # loop_id_type_pair_to_interaction[centroid_id]['basepair'][pair] = 'cWW'
            i += 1

        # help find incorrect interactions, confirm that outer interaction is set to cWW
        # counter = 0
        # for loop_id in loop_id_type_pair_to_interaction.keys():
        #     for pair, interaction in loop_id_type_pair_to_interaction[loop_id]['basepair'].items():
        #         print('loop %s pair %s apparent interaction %s' % (loop_id,pair,interaction))
        #     counter += 1
        #     if counter > 5:
        #         break

        varna_command_info, consensus_interactions = motif_to_varna(loop_id_to_core_units, loop_id_type_pair_to_interaction, cluster_method)

        print('Consensus interactions %s' % consensus_interactions)

        varna_command = 'java -cp bin/VARNAv3-7.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN "%s" -structureDBN "%s" -baseNum "#334455" -periodNum 1 -auxBPs "%s" -o %s/2ds/Group_%03d.svg' % (varna_command_info['sequenceDBN'], varna_command_info['structureDBN'], varna_command_info['auxBPs'], output_dir, group_num)
        os.system(varna_command)
        varna_command = 'java -cp bin/VARNAv3-7.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN "%s" -structureDBN "%s" -baseNum "#334455" -periodNum 1 -auxBPs "%s" -o %s/2ds/Group_%03d.png' % (varna_command_info['sequenceDBN'], varna_command_info['structureDBN'], varna_command_info['auxBPs'], output_dir, group_num)
        os.system(varna_command)
        with open(output_dir + '/' + 'varna_commands.txt', 'a+') as file:
            file.write(varna_command + '\n')

        # if VARNA is not successful, store the command and put a generic image
        varna_image_path = '%s/2ds/Group_%03d.png' % (output_dir, group_num)
        if not os.path.exists(varna_image_path):
            # write the bad command to a file so we can inspect it and fix it
            # in one case structureDBN did not have matching parentheses
            with open(output_dir + '/' + 'problemtic_varna_commands.txt', 'a+') as file:
                file.write(varna_command + '\n')
            # put a placeholder image for that group
            first_varna_image_path = '/usr/local/pipeline/hub-core/pymotifs/motif_atlas/no_varna_image_available.png'
            shutil.copy(first_varna_image_path, varna_image_path)

        varna_image_path = '%s/2ds/Group_%03d.svg' % (output_dir, group_num)
        if not os.path.exists(varna_image_path):
            # write the bad command to a file so we can inspect it and fix it
            # in one case structureDBN did not have matching parentheses
            with open(output_dir + '/' + 'problemtic_varna_commands.txt', 'a+') as file:
                file.write(varna_command + '\n')
            # put a placeholder image for that group
            first_varna_image_path = '/usr/local/pipeline/hub-core/pymotifs/motif_atlas/no_varna_image_available.svg'
            shutil.copy(first_varna_image_path, varna_image_path)

        # the following function finds and writes the bp_signature for each group using loop_id_to_core_units and consensus_interactions
        writeCSVOutput(group_num, reordered_loop_ids, loop_ids_to_discrepancy, loop_id_to_core_units, centroid_id, loop_id_type_pair_to_interaction, queries, output_dir, consensus_interactions)

        # bpSignature = find_bp_signature(queries, loop_id_to_core_units, centroid_id, loop_id_type_pair_to_interaction)
        # writeHTMLOutput(group_num,loops_of_interest, cluster_method,search_results,motif_groups,loop_id_to_core_units,reordered_loop_ids,queries,loop_ids_to_discrepancy,loop_id_to_matching_groups,loop_id_to_annotation,loop_id_to_group_number,group_max_distance,release, loop_id_type_pair_to_interaction,varna_command_info,bpSignature)

    # print(crashnow)

    print(myTimer("summary"))

