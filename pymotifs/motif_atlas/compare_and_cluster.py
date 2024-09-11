# RNAPro
# "C:\Program Files\Python38\python" compare_and_cluster.py

#from jmespath import search
from operator import truediv
from test_motif_atlas_code import add_bulged_nucleotides # this line clutters the output

from collections import defaultdict, OrderedDict
from sys import path
from sys import maxsize
from time import sleep
from time import time
import json
# import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import os.path
# import pandas as pd
import pickle
import sys
import copy
import scipy.io as sio
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

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

# this imports everything from "search.py", most notably the `function FR3D_search()`
# from search import *

from search import FR3D_search, lookUpInteractions
#from discrepancy_flip import matrix_discrepancy_cutoff_flip as matrix_discrepancy_cutoff
from fr3d_configuration import DATAPATHUNITS, DATAPATHRESULTS, DATAPATHALIGNMENTS, MOTIFPATH, DATAPATHATLAS, DATAPATHLOOPS, SERVER
from fr3d_interactions import get_fr3d_pair_to_interaction_list
from file_reading import readNAPairsFile
from file_reading import readNAPositionsFile
from file_reading import readProteinPositionsFile
from myTimer import myTimer
from orderBySimilarity import treePenalizedPathLength
# using these to get other args to feed FR3D_search()
# from pair_processing import get_pairlist # this got moved into search.py
from query_processing import calculateQueryConstraints
from query_processing import retrieveQueryInformation
from chain_to_rfam_family import read_equiv_class_csv_into_dict # this is to get quality rank of chains
from clustering_utilities import *


# from query_processing import emptyInteractionMatrix
# a program on the server can generate the next line
# That lists all the loops in a motif atlas and the unit ids in each strand.  I think.
# Keys are loop ids, then positions and maybe an indication of what strand they are in
# (1, '6DVK|1|H|G|28') means that G28 is on the border of a single-stranded region
# Then A29,
# http://rna.bgsu.edu/rna3dhub/loops/view/IL_6DVK_003


# this is where the loops are actually imported, which is not a great system
# from J3_3_70_loops_and_strands import loops_and_strands   # testing!  brand new!
# from HL_3_57_loops_and_strands import loops_and_strands
# from IL_3_57_loops_and_strands import loops_and_strands
# from IL_4_25_loops_and_strands import loops_and_strands
# the above were commented out to run on pipeline

# Chu test motifToVARNA
from motifToVARNA import get_pdb_id, find_bp_signature, get_consensus_interactions, motif_to_varna
import shutil
import csv

"""
All stacks and basepairs:
"""
bptypes = {'cWW', 'tWW', 'cWH', 'tWH', 'cWS', 'tWS', 'cHH', 'tHH', 'cHS', 'tHS',\
                   'cSS', 'tSS','cHW','tHW','cSW','tSW','cSH','tSH'}
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
"""
Download current IL json file
Save as a pickle file.
"""

def calculate_discrepancy_from_alignment(Q, ifedata, timerData):

    # this is used to avoid extra work for FR3D_search when we already have an alignment available
    # to be passed in.
    # Adam, read more about this block in FR3D_search

    # perm comes from listOfPairs from ifedata
    # perm = range(0, 48 + 1) # positions of non-bulged and aligned nucleotides
    perm = range(0, len(ifedata['centers'])) # positions of non-bulged and aligned nucleotides
    # perm is saying which nucleotides to use
    # might want a list
    # might want to pass in perm. which is a list of indices to get compared
    # that list perm would say which indices of nts from ifedata will be used
    possibility = range(0, Q['numpositions']) # set for straight alignment
    possibilities = [possibility]
    # that list possibility would say which indices of nts from Q will be used
    # might want possibility passed it

    units = ifedata['units']
    index_to_id = ifedata['index_to_id']
    pairToInteractions = ifedata['pairToInteractions']
    pairToCrossingNumber = ifedata['pairToCrossingNumber']
    # querycenters = [Q["centers"][i] for i in perm]
    # queryrotations = [Q["rotations"][i] for i in perm]
    querycenters = [Q["centers"][i] for i in possibility]
    queryrotations = [Q["rotations"][i] for i in possibility]
    timerData = myTimer("Discrepancy from query")

    possibility_to_discrepancy = {}
    for possibility in possibilities:
        possibilitycenters = []
        for i in range(0, Q['numpositions']):
            possibilitycenters.append(units[possibility[i]]["centers"])
        possibilityrotations = []
        for i in range(0, Q['numpositions']):
            possibilityrotations.append(units[possibility[i]]["rotations"])

        d = matrix_discrepancy_cutoff(querycenters, queryrotations, possibilitycenters,
            possibilityrotations, Q["discrepancy"])

        if d is not None and d < Q["discrepancy"]:
            # possibility_to_discrepancy[possibility] = d
            this_discrepancy = d

    possibilities = list(possibility_to_discrepancy.keys())

    candidates = []

    for possibility in possibilities:
        # d = possibility_to_discrepancy[possibility]
        newcandidate = {}
        indices = list(possibility)
        newcandidate['indices'] = indices
        newcandidate['unitids'] = [index_to_id[index] for index in indices]
        newcandidate['chainindices'] = [units[index]["chainindex"] for index in indices]
        newcandidate['centers'] = [units[index]["centers"] for index in indices]
        newcandidate['rotations'] = [units[index]["rotations"] for index in indices]
        # newcandidate['discrepancy'] = d
        newcandidate['discrepancy'] = this_discrepancy
        newcandidate['interactions'] = lookUpInteractions(Q,indices,
            pairToInteractions, pairToCrossingNumber, units)
        candidates.append(newcandidate)

    return(Q, candidates, timerData)


def strandify(loops_and_strands,all_structures=None):
    # loops=[]
    # for motif_group in loops_and_strands:
    #     loops_alignments = motif_group["alignment"]
    #     chainbreak = int(motif_group["chainbreak"])
    #     for loop_id,alignment in loops_alignments.items():
    #         pdb_id = loop_id.split("_")[1]
    #         if pdb_id in all_structures:
    #             new_loop = {}
    #             new_loop["loop_id"]= loop_id
    #             new_loop["strand"] = []
    #             new_loop["strand"].append(alignment[0:chainbreak])
    #             new_loop["strand"].append(alignment[chainbreak:len(alignment)])
    #             loops.append(new_loop)
    """
    The code above are for json file...
    """
    loops = []
    for loop_id in loops_and_strands:
        fields = loop_id.split("_")
        if fields[1] in all_structures:
            new_loop = {}
            new_loop["loop_id"] = loop_id
            new_loop["strand"] = [] # will be a list of 1 strand for HL, 2 for IL, 3 for J3, etc.

            current_strand = []
            bordercount = 0
            positions = sorted(loops_and_strands[loop_id].keys())
            # Loop over positions in the loop and identify where strands stop and start
            # "border" variable is 1 if the nucleotide starts or ends a single-stranded region, o/w 0
            for position in positions:
                border = loops_and_strands[loop_id][position][0]
                unit_id = loops_and_strands[loop_id][position][1]
                # print("Loop %s position %s has border %s and unit id %s" % (loop_id, position, border, unit_id))
                current_strand.append(unit_id)
                bordercount += border
                # when you get to the second bordering nucleotide on the strand, store the strand, start a new one
                if bordercount == 2:
                    bordercount = 0
                    new_loop["strand"].append(current_strand)
                    current_strand = []
                # append to the current strand until
            loops.append(new_loop)
    return loops


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


def startup_list_of_dictionaries(loops_and_strands,loops_of_interest=None):
    """
    Loop over the loop ids and extract the PDB IDs and store them in a set
    maybe also pass in "molecule type" for use in Q['requiredMoleculeType']
    """
    all_structures =[]

    if loops_of_interest:
        for loop_id in loops_of_interest:
            fields = loop_id.split("_")
            all_structures.append(fields[1])
    else:
        for loop_id in loops_and_strands.keys():
            fields = loop_id.split("_")
            all_structures.append(fields[1])

    all_structures = set(all_structures)

    pair_to_interaction_list = defaultdict(list)

    for pdb_id in all_structures:
        pair_to_interaction, pair_to_crossing = get_fr3d_pair_to_interaction_list(pdb_id)
        pair_to_interaction = strip_symmetry_from_pair_to_interaction(pair_to_interaction)

        pair_to_interaction_list.update(pair_to_interaction)

    # Loop over each loop id to identify and store the unit ids in each strand
    # loops will be a list of dictionaries, one for each loop
    # This adds a key called "strand"

    loops = strandify(loops_and_strands, all_structures)
    loops = add_bulged_nucleotides(loops, pair_to_interaction_list)

    return loops, pair_to_interaction_list

def make_query_structure(loop):
    """
    this function should take in a loop-like-object of n strands,
    and return a query dictionary `Q` (as seen in query_definitions.py)
    Here, loop is usually the strands of the **unbulged** units of the loop
    """

    positions = get_nt_positions(loop)

    Q = defaultdict(dict) # prevents later errors on assignments

    # philosophy = if the key,value takes more than one line to assign, define a `make` function for it
    Q['type'] = "mixed"     # geometric AND symbolic
    Q['name'] = "AVA"       # this is to DODGE Q possibly getting re-defined # could be named by loop ID

    Q['motif_atlas'] = True # when this field exists, FR3D search will not print list lengths

    if len(loop['strand']) <= 2:
        Q['discrepancy'] = 1.0  # 1.0 is what the Motif Atlas uses for HL and IL
    else:
        Q['discrepancy'] = 1.0  # in case we want to try something new for J3

    Q = make_num_positions(Q, positions)

    # seeding section
    Q["interactionMatrix"] = emptyInteractionMatrix(Q['numpositions'])
    Q['unitID'] = []        # these are the unit ids that the query will actually use for the search
    Q['requiredUnitType'] = [None] * Q['numpositions'] # these list operations are copied from query_processing.py
    Q['requiredMoleculeType'] = ["RNA"] * Q['numpositions'] # future, find a way to discern in DNA/RNA, or be a user input variable

    Q['errorMessage'] = []
    Q['MAXTIME'] = 60*60    # maximum search time 1 hour, may not be imposed
    Q["CPUTimeUsed"] = 0    # needs to be tracked, but actually ignored here

    # below are things that prefer to be fed strand by strand
    for index, strand in enumerate(loop['strand']):
        Q = make_unit_ids(Q, strand)
        if len(positions[index]) > 1:
            # if strand is more than 1 nt, find direction + constraint, else `pass`
            Q = make_directional_constraints(Q, positions[index], "increasing")

    # below are things requiring above functions to be done, but dislike strand by strand
    Q = make_search_files(Q)

    if len(positions) == 0:
        print("344 make_query_structure: No nt positions in loop %s" % loop)
        Q['errorStatus'] == "write and exit"
        return Q
    elif len(positions[0]) == 0:
        print("348 make_query_structure: No nt positions first strand of loop %s" % loop)
        Q['errorStatus'] == "write and exit"
        return Q

    Q = make_req_interactions(Q, positions)

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

    if len(positions) == 0:
        print("370: No nt positions in loop %s" % loop)
        Q['errorStatus'] == "write and exit"
        return Q
    elif len(positions[0]) == 0:
        print("374: No nt positions first strand of loop %s" % loop)
        Q['errorStatus'] == "write and exit"
        return Q

    for i in range(len(positions)):
        flanking_positions[i] = [ positions[i][0],positions[i][-1]]

    # we probably need a new function here to return only the flanking nt positions
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

    Q["interactionMatrix"] = emptyInteractionMatrix(Q['numpositions'])

    #Only for IL
    if Q['numpositions'] == 4:
        Q["interactionMatrix"][0][3] = "cWW and GC CG AU UA GU UG"
        Q["interactionMatrix"][1][2] = "cWW and GC CG AU UA GU UG"
        Q["interactionMatrix"][1][0] = ">"
        Q["interactionMatrix"][3][2] = ">"

    if Q['numpositions'] == 6:
        Q["interactionMatrix"][0][5] = "cWW and GC CG AU UA GU UG"
        Q["interactionMatrix"][1][2] = "cWW and GC CG AU UA GU UG"
        Q["interactionMatrix"][3][4] = "cWW and GC CG AU UA GU UG"
        Q["interactionMatrix"][1][0] = ">"
        Q["interactionMatrix"][3][2] = ">"
        Q["interactionMatrix"][5][4] = ">"

    Q['unitID'] = []
    Q['requiredUnitType'] = [None] * Q['numpositions'] # these list operations are copied from query_processing.py
    Q['requiredMoleculeType'] = ["RNA"] * Q['numpositions'] # future, find a way to discern in DNA/RNA, or be a user input variable

    Q['errorMessage'] = []

    # below are things that prefer to be fed strand by strand
    for index, strand in enumerate(flanking_strands):
        Q = make_unit_ids(Q, strand)
        # if len(flanking_positions[index]) > 1:
        # # if strand is more than 1 nt, find direction + constraint, else `pass`
        #     Q = make_directional_constraints(Q, flanking_positions[index], "increasing")

    # below are things requiring above functions to be done, but dislike strand by strand
    Q = make_search_files(Q)

    # this will have fewer interactions when just searching for flanking pairs
    #Q = make_req_interactions(Q, positions)

    return Q


def get_nt_positions(loop):
    """
    takes in a list of strings representing nucleotides
    returns a list of lists of nucleotide position numbers (by strand)
    tiny ex: positions = [[nt1, nt2, nt3], [nt4, nt5, nt6]]
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

    return(positions)


def make_directional_constraints(query, positions, order):
    """
    Adam, make a good description of this
    """

    if order.lower() == "increasing": # these LOOK backwards, but are not
        symbol = ">"
    elif order.lower() == "decreasing":
        symbol = "<"
    else:
        raise NameError('Tried to make directional without "increasing" or "decreasing" order argument')

    previousPositions = [] # needed to FILL lower diagonal
    firstIter = True
    for position in positions: # remember that this function is called on each strand
        if firstIter == False: # if first iteration, previousPositions is empty
            for prev in previousPositions:
                query['interactionMatrix'][position][prev] = symbol
            previousPositions.append(position)
        else:
            firstIter = False
            previousPositions.append(position)

    return(query)


def make_num_positions(query, positions):
    """
    counts the nts in a loop
    """

    # adam, consider removing this
    counter = 0
    for strand in positions:
        counter += len(strand)
    query['numpositions'] = counter

    return(query)

def emptyInteractionMatrix(n):
    """set up an empty list of required interactions,
    then the user only needs to set the actual constraints needed"""

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
            query['interactionMatrix'][prevTail][head] = "cWW GC CG AU UA GU UG"
        if tail == tails[-1]: # if last pass
            query['interactionMatrix'][firstHead][tail] = "cWW GC CG AU UA GU UG"
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

    return(query)


def make_unit_ids(query, strand):
    """example of what this should make
    Q["unitID"] = ["4V9F|1|0|U|1026","4V9F|1|0|A|1032","4V9F|1|0|G|1034", ... ]"""

    for nt in strand:
        query['unitID'].append(nt)

    return(query)


def combine_dicts(x, y):
    z = x.copy()
    z.update(y)
    return(z)


def readPAI(query, ifename, alternate = "", by_chain = False):
    """
    WAS readPositionsAndInteractions() from ifedata.py.
    It got copied + changed here for AvA needs
    """

    RNA_positions_file_name = ifename.replace("|","-").replace("+","_")
    starting_index = 0
    ifedata = {}
    ifedata['index_to_id'] = {}
    ifedata['id_to_index'] = {}
    ifedata['centers'] = np.empty((0, 3))
    ifedata['rotations'] = np.empty((0, 3, 3))
    ifedata['ids'] = []
    ifedata['units'] = []
    allCenters = []
    allModels = []

    # lists should start empty, append with RNA if necessary, with protein if necessary, with DNA if necessary, etc.
    # if(any("RNA" in query["requiredMoleculeType"][index] for index in range(len(query["requiredMoleculeType"])))):
        # check to see if RNA is a required unit type, and if so, read RNA data

    for chain_string in RNA_positions_file_name.split("_"):
        # print("in readPAI reading %s" % (chain_string))
        query, centers, rotations, ids, id_to_index, index_to_id, chainIndices = readNAPositionsFile(query, chain_string, starting_index)

        if len(centers) > 0:
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
                allModels.append(data[1]) # extract model number
                unit_information["unitType"] = data[3] # supposed to get the nucleotide letter
                unit_information["moleculeType"] = "RNA"
                unit_information["chainindex"] = chainIndices[id_to_index[unitID]-starting_index]
                ifedata["units"].append(unit_information)

            starting_index += len(centers)


    PDBID = ifename.split("|")[0]
    query, interactionToPairs, pairToInteractions, pairToCrossingNumber = readNAPairsFile(query, PDBID, ifedata["id_to_index"], alternate)
    ifedata['interactionToPairs'] = interactionToPairs
    ifedata['pairToInteractions'] = pairToInteractions
    ifedata['pairToCrossingNumber'] = pairToCrossingNumber
    ifedata['models'] = allModels

    return query, ifedata

# Read .pickle file of RNA base center and rotation matrix; download if necessary
def read_RNA_pos_file(Q, chain_string, starting_index, by_chain = False):

    ids = []
    chainIndices = []
    centers = []
    rotations = []
    line_num = starting_index
    id_to_index = defaultdict()
    index_to_id = defaultdict()
    our_indices = []

    filename = chain_string + "_RNA.pickle"
    filename = chain_string + "_NA.pickle"   # calculate discrepancy using glycosidic atom

    pathAndFileName = os.path.join(DATAPATHUNITS, filename)
    print("pathandfilename: %s" % (pathAndFileName))

    if not os.path.exists(DATAPATHUNITS):
        os.mkdir(DATAPATHUNITS)

    if not os.path.exists(pathAndFileName) and not SERVER:
        print("Downloading "+filename)
        urlretrieve("http://rna.bgsu.edu/units/" + filename, pathAndFileName)


    if os.path.exists(pathAndFileName):

        print("read_RNA_pos_file: Reading "+pathAndFileName)

        if sys.version_info[0] < 3:
            try:
                ids, chainIndices, centers, rotations = pickle.load(open(pathAndFileName,"rb"))
            except:
                print("Could not read "+filename+" A==================================")
                Q["userMessage"].append("Could not retrieve RNA unit file "+filename)
                centers = np.zeros((0, 3))  #of lines/nucleotides in file
                rotations = np.zeros((0, 3, 3))
        else:
            try:
                ids, chainIndices, centers, rotations = pickle.load(open(pathAndFileName,"rb"), encoding = 'latin1')
            except:
                print("Could not read "+filename+" A*********************************")
                Q["userMessage"].append("Could not retrieve RNA unit file "+filename)
                centers = np.zeros((0, 3))  # of lines/nucleotides in file
                rotations = np.zeros((0, 3, 3))


        if(by_chain == False):
            for i in range(0, len(ids)):
                if ids[i] in Q['fullUnits']:
                    our_indices.append(i) # above line could also add "or in Q['bulges']"
                    id_to_index[ids[i]] = line_num
                    index_to_id[line_num] = ids[i]
                    line_num += 1
        else: # if just going by a WHOLE CHAIN at a time
            for i in range(0, len(ids)):
                our_indices.append(i) # above line could also add "or in Q['bulges']"
                id_to_index[ids[i]] = line_num
                index_to_id[line_num] = ids[i]
                line_num += 1

    else:
        print("Could not find "+filename)
        Q["userMessage"].append("Could not retrieve RNA unit file "+filename)
        centers = np.zeros((0, 3))  # of lines/nucleotides in file
        rotations = np.zeros((0, 3, 3))

    ids = [ids[i] for i in our_indices] # equivalent to: ids = ids[our_indices]
    centers = np.asarray([centers[i] for i in our_indices])
    rotations = np.asarray([rotations[i] for i in our_indices])
    chainIndices = [chainIndices[i] for i in our_indices]

    return Q, centers, rotations, ids, id_to_index, index_to_id, chainIndices


def create_query(loop, ifedata = {}):
    """
    Make a query data structure based on the data in loop
    Makes the query (Q) without bulges, while bulges should be left IN targets (ifedata)
    """

    query_length = sum([len(strand) for strand in loop["strand"]])
    query_bulged_length = len(loop["bulged"])

    # identify non-bulged nucleotides to form the basis of the query
    # double list comprehension to keep the structure of ['strand'] being a list of lists, while removing bulges
    unbulged = {}
    unbulged['strand'] = [[unitID for unitID in strand if unitID not in loop['bulged']] for strand in loop['strand']]

    # when a 5-nt IL with one bulged nt, keep all 5 nucleotides in unbulged variable
    if query_length == 5 and query_bulged_length == 1:
        query = make_query_structure(loop)
        query['discrepancy']= 98
    else:
        query = make_query_structure(unbulged)

    if 'errorStatus' in query and query['errorStatus'] == "write and exit":
        return None

    query['fullUnits'] = [unitID for strand in loop['strand'] for unitID in strand]
    if not ifedata:
        print("No ifedata in %s" % loop["loop_id"])
        return None
    else:
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

    return query


def create_search_space(loop):
    query = defaultdict(list)
    # makes the query (Q) without bulges, while bulges should be left IN targets (ifedata)
    for strand in loop['strand']:
        for nt in strand:
            query['unitID'].append(nt)
    positions = get_nt_positions(loop)
    query = make_num_positions(query, positions)

    query['fullUnits'] = [unitID for strand in loop['strand'] for unitID in strand]
    query['requiredMoleculeType'] = ["RNA"] * query['numpositions']
    query['activeInteractions'] = ['cWW']
    query['userMessage']=[]
    query = make_search_files(query)

    print()
    print("inside create search space")
    print(query)
    print(query['searchFiles'][0])
    print("entering read PAI")
    query, ifedata = readPAI(query, query['searchFiles'][0])
    print()
    print("end create search space")
    print()
    return(ifedata)



def create_flanking_bp_query(loop, ifedata = {}):
    unbulged = {}

    unbulged['strand'] = [[unitID for unitID in strand if unitID not in loop['bulged']] for strand in loop['strand']]
    flanking_bp_query = make_flanking_bp_query_structure(unbulged)

    if not flanking_bp_query:
        return None

    flanking_bp_query = retrieveQueryInformation(flanking_bp_query, ifedata)

    if 'errorStatus' in flanking_bp_query and flanking_bp_query['errorStatus'] == "write and exit":
        return None
    elif not flanking_bp_query:
        return None
    else:
        flanking_bp_query = calculateQueryConstraints(flanking_bp_query)
        flanking_bp_query['numStrands'] = len(loop['strand'])
    return(flanking_bp_query)


def get_bulge(query,search_space):
    query_bulge = ''
    search_space_bulge = ''
    if len(query['loop_info']['bulged'])>0:
        query_bulge = query['loop_info']['bulged'][0].split("|")[3]
    if len(search_space['loop_info']['bulged'])>0:
        search_space_bulge = search_space['loop_info']['bulged'][0].split("|")[3]
    return(query_bulge,search_space_bulge)


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


def filter_out_conflicting_basepairs_and_stacks(candidates, query, search_space):
    """
    Look for conflicting base pairs or stacking between query and candidates
    Return candidates with no conflicts
    candidates data is just the ifedata of the candidate
    """

    query_ifedata = query['ifedata']

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
            index_i_of_query = query_ifedata['id_to_index'][ query["Q"]["unitID"][i] ]
            index_i_of_candidate = candidate['indices'][i] #The integer indentifying nucleotide in i position
            matched_search_space_indices.add(index_i_of_candidate)

            for j in range(i + 1, len(candidate['indices'])): # goes like i+1 to 7
                index_j_of_query = query_ifedata['id_to_index'][ query["Q"]["unitID"][j] ]
                index_j_of_candidate = candidate['indices'][j]

                query_interaction = get_set_of_interactions(index_i_of_query, index_j_of_query, query_ifedata) & (bptypes | stacks)
                candidate_interaction = get_set_of_interactions(index_i_of_candidate, index_j_of_candidate, candidates_data) & (bptypes | stacks)
                DQ_code += are_conflicting(query_interaction, candidate_interaction)

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

    return not_rejected


def keep_lowest_discrepancy_candidate(candidates):
    """
    this takes a list of possible candidates and returns the one
    with the lowest discrepancy, unless all candidates have
    disqualification codes, then it returns an empty list
    """
    lowest_dq_discrepancy = 100 # purposely seeded high
    lowest_matched_discrepancy = 100
    positive_found = 0
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
    if((i,j) in ifedata['pairToInteractions'] ):
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

    if(q_int == c_int):
        return []
    if( (q_int & stacks) and (c_int & stacks) ): # loose due to modeling flips
        return []
    if( not ( q_int&stacks) and (c_int&stacks) ):
        if(len(q_int) != 0):
            return [BASESTACK_MISMATCH] # base stack mismatch
        else:
            return []
    if( (q_int & stacks) and not (c_int & stacks)):
        if(len(c_int) != 0):
            return [BASESTACK_MISMATCH] # base stack mismatch
        else:
            return []
    # no stacks by this line
    if( len(q_int) == 0 or len(c_int) == 0 ): # one is empty
        return [] # compatible
    if( q_int != c_int ): # i dont believe any loops have made it this far yet (82x82 tested)
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
        return(results)
    #file_name_and_path = path + "/search_results/" + structure_name + ".pickle"
    file_name = loop_type + "_" + pdb_id + "_search_results.pickle"
    file_path = os.path.join(path,file_name)

    if os.path.exists(file_path):
        if sys.version_info[0] < 3:
            results = pickle.load(open(file_path, "rb"))
        else:
            results = pickle.load(open(file_path, "rb"), encoding = 'latin1')
        # print("debug info: load previous search results: %s" %results)
        for (query,search_space),value in results.items(): #Load No-match result as [{'dq': 0, 'discrepancy': 99}] instead of an integer
            if isinstance(value,int):
                results[(query,search_space)] = [{'dq': [value], 'discrepancy': 99}]

    return(results)


def save_search_results_one_pdb(loop_type, pdb_id, results, path = DATAPATHRESULTS):
    if not os.path.exists(path):
        os.mkdir(path)

    file_name = loop_type + '_' + pdb_id + '_search_results.pickle'
    file_path = os.path.join(path, file_name)

    pickle.dump(obj = results, file = open(file_path, "wb"), protocol = 2)
    return()

###### This is unused, just a stub
def look_for_a_structure_alignment(structure_a, structure_b):
    """
    in the future, this will likely call an API or use a
    function from Satya
    """
    found = False
    alignment = {"loop to": "loop"}
    if found:
        return(alignment)
    else:
        return(False)


def list_all_chains_in(loops_from_structure, best_quality_chains = False):
    from itertools import permutations
    accumulator_of_chains = []

    # this top "dict" block... I'm not sure if it's used anymore (09/18/2023 - Adam)
    if type(loops_from_structure) is dict:
        for loop_id, loop_info in loops_from_structure.items():
            this_loop_chains = []
            this_loop_nt_ids = [nt for strand in loop_info['loop_info']['strand'] for nt in strand]
            for nt_id in this_loop_nt_ids:
                fields = nt_id.split("|")
                one_chain = "|".join(fields[0:3])
                this_loop_chains.append(one_chain)
            accumulator_of_chains = accumulator_of_chains + ["+".join(set(sorted(this_loop_chains)))]


    # if "loop_info" in loops_from_structure:
    if type(loops_from_structure) is list:
        for loop_dict in loops_from_structure:
            this_loop_chains = []
            this_loop_nt_ids = [nt for strand in loop_dict['strand'] for nt in strand]
            for nt_id in this_loop_nt_ids:
                fields = nt_id.split("|")
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
    if not os.path.exists(DATAPATHALIGNMENTS):
        os.mkdir(DATAPATHALIGNMENTS)

    # if we have the alignment. load it and move on
    if(os.path.exists(path_and_file)):
        ntB_to_ntA_alignments = pickle.load(open(path_and_file, 'rb'))
        return(ntB_to_ntA_alignments)

    print("try to open: " + "http://rna.bgsu.edu/correspondence/align_chains?chains=" +
        chain_of_b + "," + chain_of_a)
    try:
        reply = urlopen("http://rna.bgsu.edu/correspondence/align_chains?chains=" +
        chain_of_b + "," + chain_of_a).read()
        print("success")
        reply = reply.decode("ascii")
    except:
        reply = ""
        print("failed")

    for index, line in enumerate(reply.split("\n")):
        fields = line.split()
        if(len(fields) == 2):
            ntB, ntA = fields[0], fields[1]
            if ntB == "NULL" or ntA == "NULL":
                pass
            else:
                # there ARE possible overwrites here =/
                ntB_to_ntA_alignments[ntB] = ntA


    # save this as a pickle file for next time
    with open(path_and_file, 'wb') as f:
        pickle.dump(ntB_to_ntA_alignments, f)

    return(ntB_to_ntA_alignments)


def map_the_loops_of_these_chains(queries, search_spaces, q_chains, ss_chains, nt_to_loop_of_ss):
    """
    this is mapping loops to loops from two different chains in the same equivalence class
    """
    # `ss` is short for search space, `q` is short for query

    mapping = {}
    ntB_to_ntA_alignment = {}

    # ntB_to_ntA_alignments[structure a nt unit id] = structure b nt unit id

    for q_chain in q_chains.split("+"):
        for ss_chain in ss_chains.split("+"):
            ntB_to_ntA_alignment.update(get_nt_alignment_between_chains(chain_of_a = ss_chain, chain_of_b = q_chain))

    if len(ntB_to_ntA_alignment) == 0:
        return()

    # ^ could be .update() instead to be passed in and returned, also needs chain_to_equiv_class dict


    for loop_name_in_q, q_loop_data in queries.items():
        # query_loop_data['loop_info']['strand']] <- this is a list of lists, ex: [["6DVK|1|H|G|28", "6DVK|1|H|A|29"], ["6DVK|1|H|G|61", "6DVK|1|H|C|62"]]
        # when a nt from search_space is aligned to a nt in a query loop, that loop's name is appended

        # this next line is a list comprehension, to flatten a list of lists into a regular 1D list
        nts_of_q_loop = [nt for strand in q_loop_data['loop_info']['strand'] for nt in strand]
        ss_loops_with_aligned_nts = []

        for nt in nts_of_q_loop:
            # new
            ss_loop = nt_to_loop_of_ss.get(ntB_to_ntA_alignment.get(nt))
            if(ss_loop):
                ss_loops_with_aligned_nts.append(ss_loop)
            # old
            # if nt in ntB_to_ntA_alignment:
            #   if ntB_to_ntA_alignment[nt] in nt_to_loop_of_ss:
            #       ss_loops_with_aligned_nts.append(nt_to_loop_of_ss[ntB_to_ntA_alignment[nt]])

        # next line find the most common match
        # max(set(query_loop_matches), key = query_loop_matches.count)
        if(ss_loops_with_aligned_nts):
            most_common_ss_loop = max(set(ss_loops_with_aligned_nts), key = ss_loops_with_aligned_nts.count)

            # is the most common match over HALF of the nts?
            if(ss_loops_with_aligned_nts.count(most_common_ss_loop)/len(ss_loops_with_aligned_nts) > .5):
                # add it to the loop to loop mapping for return
                mapping[loop_name_in_q] = most_common_ss_loop


    return(mapping)


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


def run_aligned_searches_first(queries = {}, search_spaces = {}, structure_a_best_chains = False, flanking_bp_queries = {}, return_all_queries = False, force_alignment = False):
    """
    this is calling fr3d on a loop-by-loop basis
    where matches are almost guaranteed
    """
    # this function still needs work Adam
    # queries = structure B, search_spaces = structure A

    # aligned_search_pairs[query_loop_id] = search_space_loop_id
    aligned_search_cutoff = .6
    aligned_search_pairs = {}
    # queries_to_still_run = {}
    search_spaces_without_matches = {}
    search_results = {}
    failed_aligned_search_pairs = {}
    Q = {}

    # could be a way to focus down list of chains in a
    query_chains = list_all_chains_in(queries)

    # best chains stuff should be done before this function is called
    if(structure_a_best_chains):
        search_space_chains = structure_a_best_chains
    else:
        search_space_chains = list_all_chains_in(search_spaces)

    # chain_to_equiv_class
    # chain_to_equiv_class = read_equiv_class_csv_into_dict()
    # breakpoint()
    nt_to_loop_of_ss = nt_to_loop_mapping(search_spaces)
    nt_to_loop_of_q = nt_to_loop_mapping(queries)

    # filter chain list of search space to "best chains"
    # some_new_function_to_do_that()

    if not os.path.exists(DATAPATHALIGNMENTS):
        os.mkdir(DATAPATHALIGNMENTS)
    files_of_dir = os.listdir(DATAPATHALIGNMENTS)
    txt_files = [file for file in files_of_dir if ".txt" in file]

    for ss_chain in search_space_chains:
        for q_chain in query_chains:

            # if force_alignment:
            # map the loops of these chains will need to return some alignment data
            priority_searches = map_the_loops_of_these_chains(queries, search_spaces, q_chain, ss_chain, nt_to_loop_of_ss)
            aligned_search_pairs.update(priority_searches)
            # else: # without force alignment, check if equivalence classes match
            #     if(chain_to_equiv_class.get(q_chain)):
            #         if(chain_to_equiv_class.get(ss_chain)):
            #             if(chain_to_equiv_class[q_chain]["equivalence class"] == chain_to_equiv_class[ss_chain]["equivalence class"]):
            #                 priority_searches = map_the_loops_of_these_chains(queries, search_spaces, q_chain, ss_chain, nt_to_loop_of_ss)
            #                 aligned_search_pairs.update(priority_searches)
    timer_data = myTimer("Running Aligned Searches")


    for query_id, search_space_id in sorted(aligned_search_pairs.items()):
        print("Aligned search %s in %s" % (query_id, search_space_id))
        #print("Query:",queries[query_id]['Q']['unitID'])
        queries[query_id]['Q']['discrepancy'] = aligned_search_cutoff

        # for i in range(0, queries[query_id]['Q']['numpositions']):
        #   queries[query_id]['Q']['cutoff'][i] = aligned_search_cutoff

        search_space_length = len(search_spaces[search_space_id]['ifedata']['centers'])
        candidates = []


        # if we can align directly from query (or query with bulges) to search space, use this shortcut
        if queries[query_id]['Q']['numpositions'] == search_space_length:
            Q, candidates, elapsed_time = calculate_discrepancy_from_alignment(Q = queries[query_id]['Q'],
                ifedata = search_spaces[search_space_id]['ifedata'],
                timerData = timer_data)
        elif len(queries[query_id]['Q']['fullUnits']) == search_space_length:
            # if the bulges are important
            mockQ = copy.deepcopy(queries[query_id]['Q'])
            mockQ['numpositions'] = len(mockQ['fullUnits'])
            mockQ['centers'] = queries[query_id]['ifedata']['centers']
            mockQ['rotations'] = queries[query_id]['ifedata']['rotations']
            Q, candidates, elapsed_time = calculate_discrepancy_from_alignment(Q = mockQ,
                ifedata = search_spaces[search_space_id]['ifedata'],
                timerData = timer_data)
        if not candidates: # if we skipped shortcut or shortcut was bad, do a regular fr3d search
            if search_space_length < 20:
                Q, candidates, elapsed_time = FR3D_search(Q = queries[query_id]['Q'],
                                                ifedata = search_spaces[search_space_id]['ifedata'],
                                                ifename = search_space_id,
                                                timerData = timer_data)

        if candidates: # this is to filter out empties AND get rid of all being list[0]
            temp_result = filter_out_conflicting_basepairs_and_stacks(candidates, queries[query_id], search_spaces[search_space_id])
            if temp_result:
                temp_result = filter_on_unmatched_nucleotides(temp_result,search_spaces[search_space_id],queries[query_id])
                temp_result = keep_lowest_discrepancy_candidate(temp_result)
                # add cutoff here?

                if temp_result:
                    if temp_result[0]["discrepancy"] <= aligned_search_cutoff:
                        temp_result[0]['match_type'] = "homologous"
                        search_results[(query_id,search_space_id)] = temp_result
                    else: # if discrepancy too high
                        failed_aligned_search_pairs[query_id] = aligned_search_pairs[query_id]
                else:
                    failed_aligned_search_pairs[query_id] = aligned_search_pairs[query_id]
            else: # if NOT temp_result
                failed_aligned_search_pairs[query_id] = aligned_search_pairs[query_id]
        else: # if no candidates from fr3d_search
            failed_aligned_search_pairs[query_id] = aligned_search_pairs[query_id]

    # if not return_all_queries:
    #     # queries dict minus aligned_search_pairs keys
    #     for loop_id, query_info in queries.items():
    #         if loop_id not in aligned_search_pairs:
    #             queries_to_still_run[loop_id] = query_info
    #         if loop_id in failed_aligned_search_pairs:
    #             queries_to_still_run[loop_id] = query_info
    # else:
    #     queries_to_still_run = queries

    # ss dict minus aligned_search_pairs keys
    for loop_id, ss_info in search_spaces.items():
        if loop_id not in aligned_search_pairs.values():
            search_spaces_without_matches[loop_id] = ss_info
        if loop_id in failed_aligned_search_pairs.values():
            search_spaces_without_matches[loop_id] = ss_info

    # breakpoint()
    return(search_results, search_spaces_without_matches)


# test_loop = [{'loop_id': 'IL_5J7L_001', 'strand': [['5J7L|1|AA|U|30', '5J7L|1|AA|G|31', '5J7L|1|AA|A|32'], ['5J7L|1|AA|U|552', '5J7L|1|AA|A|553']], 'bulged': ['5J7L|1|AA|G|31']}]
# test_double_chain_loop = [{'loop_id': 'IL_5TBW_002', 'strand' : [['5TBW|1|1|C|7', '5TBW|1|1|C|8'], ['5TBW|1|4|G|150', '5TBW|1|4|C|151', '5TBW|1|4|G|152']], 'bulged': ['5J7L|1|AA|G|31']}]


def create_queries_and_search_spaces(loops, loops_of_interest = None):
    queries = {}
    flanking_bp_queries = {}
    search_spaces = {}
    current_pdb_structure = ''
    pdb_data = {}
    # member_to_equiv_class = read_equiv_class_csv_into_dict()
    # chain_to_rfam = read_chain_to_rfam()
    # make sure loops is sorted by PDB id and if possible by chain
    # loops ARE sorted by loop number which seems to follow the above convention

    # sort by loop_id here anyway (lambda function)

    for loop in sorted(loops, key=lambda k: k['loop_id']):
        # loop['rfam'] = which_rfam_and_equiv(loop, chain_to_rfam, member_to_equiv_class) # keep this guy from deletion
        pdb = loop['loop_id'].split("_")[1]

        # print('Creating query for loop: %s' % loop['loop_id'])

        # print('1368 create_queries_and_search_spaces: Current loop:')
        # print(loop)

        # if we're on a new pdb, read all pairs
        if pdb != current_pdb_structure:
            pdb_data = {}
            current_pdb_structure = pdb

            # read all pairwise interactions in this structure
            pair_to_interaction_list, pair_to_crossing_number = get_fr3d_pair_to_interaction_list(pdb, near = True)

            pair_to_interaction_list = strip_symmetry_from_pair_to_interaction(pair_to_interaction_list)

        # get a flat list of the nucleotides (full unit ids) of this loop
        this_loop_nt_ids = [nt for strand in loop['strand'] for nt in strand]
        # print("this_loop_nt_ids: %s" % (this_loop_nt_ids))

        # here is test load block
        # dodge using rfam
        # no rfam loop id ['']
        this_loop_chains = []
        for nt_id in this_loop_nt_ids:
            fields = nt_id.split("|")
            one_chain = "|".join(fields[0:3])
            this_loop_chains.append(one_chain)

        chain = "+".join(set(sorted(this_loop_chains)))
        # the above line sometimes leaves a leading "+"
        # subtract below
        if chain:
            if chain[0] == "+":
                chain = chain[1:]

        if not chain in pdb_data:
            empty_query = defaultdict(list)
            empty_query, pdb_data[chain] = readPAI(empty_query, chain, by_chain = True)

        ifedata = {}
        ifedata['index_to_id'] = {}
        ifedata['id_to_index'] = {}
        ifedata['centers'] = np.empty((0, 3))
        ifedata['rotations'] = np.empty((0, 3, 3))    # is this actually needed?
        ifedata['ids'] = this_loop_nt_ids
        ifedata['units'] = []                         # heavily used by FR3D search
        ifedata['interactionToPairs'] = {}            # key will be interaction type
        ifedata['pairToInteractions'] = defaultdict(list)
        ifedata['pairToCrossingNumber'] = {}
        ifedata['models'] = []

        # convert chain_indices from the full pdb_data
        # to shallow_indices (loop indices), and extract only data
        # relevant to this loop
        for shallow_index, nt_id in enumerate(this_loop_nt_ids):
            # fields = nt_id.split("|") # these got commented out with the new load style
            # chain = "|".join(fields[0:3])

            chain_index = None  # occasionally there is no match

            if nt_id in pdb_data[chain]['id_to_index']:
                chain_index = pdb_data[chain]['id_to_index'][nt_id]
            else:
                # Sometimes unit id like 1NUV|1|A|G|26 needs to match 1NUV|1|A|G|26||A
                for check_nt_id, index in sorted(pdb_data[chain]['id_to_index'].items()):
                    if nt_id in check_nt_id:
                        remainder = check_nt_id.replace(nt_id,"")
                        if "A" in remainder:
                            chain_index = index
                            break
                        elif "B" in remainder:
                            chain_index = index
                            break
                        else:
                            chain_index = index
                            break

            if chain_index == None:
                print("create_queries_and_search_spaces: Could not find chain index for %s ==============================" % nt_id)
                print(sorted(pdb_data[chain]['id_to_index'].items()))
                ifedata['missing_center'] = True
            else:
                ifedata['id_to_index'][nt_id] = shallow_index
                ifedata['index_to_id'][shallow_index] = nt_id
                ifedata['centers'] = np.append(ifedata['centers'], np.asarray([pdb_data[chain]['centers'][chain_index]]), axis = 0)
                # consider using dummy data below to test of "rotations" is needed, Adam
                ifedata['rotations'] = np.append(ifedata['rotations'], np.asarray([pdb_data[chain]['rotations'][chain_index]]), axis = 0)
                # ifedata2['rotations'] = np.append(ifedata2['rotations'], np.asarray([[1,2,3],[4,5,6],[7,8,9]]), axis = 0)
                ifedata['units'].append(pdb_data[chain]['units'][chain_index])
                ifedata['models'].append(pdb_data[chain]['models'][chain_index])


        # new code to get the required interactions from interactionToTriple
        for shallow_index_1, nt_id_1 in enumerate(this_loop_nt_ids):
            # breakpoint()

            for shallow_index_2, nt_id_2 in enumerate(this_loop_nt_ids):
                if (nt_id_1, nt_id_2) in pair_to_interaction_list:
                    # store the interaction according to the shallow index, not just the unit ids
                    index_pair = (shallow_index_1, shallow_index_2)
                    # changed 03/15/23
                    ifedata['pairToInteractions'][index_pair] = pair_to_interaction_list[(nt_id_1,nt_id_2)]
                    ifedata['pairToCrossingNumber'][index_pair] = pair_to_crossing_number[(nt_id_1,nt_id_2)]

                    crossing_number = pair_to_crossing_number[(nt_id_1,nt_id_2)]

                    # needed by FR3D search for sure
                    for interaction in pair_to_interaction_list[(nt_id_1, nt_id_2)]:
                        if interaction in all_bptypes | all_stacks:
                            if not interaction in ifedata['interactionToPairs']:
                                ifedata['interactionToPairs'][interaction] = [[],[]]
                            ifedata['interactionToPairs'][interaction][0].append(index_pair)
                            ifedata['interactionToPairs'][interaction][1].append(crossing_number)

        # Above currently includes non all_stacks | all_bp_types
        # Debugging query_length>2:
        if not 'missing_center' in ifedata:

            # print("1494 create_queries_and_search_spaces loop")
            # print(loop)

            Q = create_query(loop, ifedata)

            # print("1499 create_queries_and_search_spaces loop")
            # print(loop)

            if not Q:
                continue
            elif "errorStatus" in Q and Q["errorStatus"] == "write and exit":
                continue

            # print("1504 create_queries_and_search_spaces Q")
            # print(Q)

            if len(Q['fullUnits']) > 2:
                bp_Q = create_flanking_bp_query(loop, ifedata)

                if not bp_Q:
                    continue

                search_spaces[loop['loop_id']] = {}
                search_spaces[loop['loop_id']]['loop_info'] = loop
                search_spaces[loop['loop_id']]['ifedata'] = ifedata

                queries[loop['loop_id']] = {}
                queries[loop['loop_id']]['loop_info'] = loop
                queries[loop['loop_id']]['Q'] = Q
                queries[loop['loop_id']]['ifedata'] = ifedata

                flanking_bp_queries[loop['loop_id']] = {}
                flanking_bp_queries[loop['loop_id']]['loop_info'] = loop
                flanking_bp_queries[loop['loop_id']]['Q'] = bp_Q
                flanking_bp_queries[loop['loop_id']]['ifedata'] = ifedata

    return queries, search_spaces, flanking_bp_queries



def create_flanking_bp_queries(loops, loops_of_interest = None, search_spaces = {}):
    print("####################################")
    print("This function SHOULD be defunct. Deleting soon")
    print("####################################")
    queries = {}
    for loop in loops:
        if not loops_of_interest or loop['loop_id'] in loops_of_interest:
            if search_spaces[loop['loop_id']]:
                Q = create_flanking_bp_query(loop, search_spaces[loop['loop_id']]['ifedata'])
            else:
                Q = create_flanking_bp_query(loop)
                print("This function 'create_flanking_bp_queries' should NOT be called without supplying search_spaces (aka ifedata)")
            if Q:
                queries[loop['loop_id']] = {}
                queries[loop['loop_id']]['loop_info'] = loop
                queries[loop['loop_id']]['Q'] = Q
                if loop['loop_id'] == "IL_5TBW_144" or loop['loop_id'] == "IL_6SY6_001":
                    queries[loop['loop_id']]['Q']['discrepancy']= 0.95
    return(queries)


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


def set_loop_type(queries):
    loop_id = list(queries.keys())[0]
    loop_type = loop_id.split('_')[0]
    return loop_type


def all_against_all_searches(loop_type,queries, search_spaces,flanking_bp_queries,load_saved_searches = True,reversed_search=False, save_path = DATAPATHRESULTS):
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
    query_pdb_id = ""
    search_space_pdb_id = ""
    search_results = {} # search_results store the final, accumulated results

    loop_counter = 0

    for search_space_id in sorted(search_spaces.keys(), reverse = reversed_search):
        loop_counter += 1
        print("Searching for %d query loops inside loop %s, %d/%d" % (len(queries),search_space_id,loop_counter,len(search_spaces)))

        # Only load previous result when there is a change in search_space_pdb_id.
        # Save the results if not the first iteration.
        if search_space_pdb_id != search_space_id.split("_")[1]:

            # moved saving up here to only happen once per pdb
            if search_space_pdb_id != "":
                if new_results:
                    save_search_results_one_pdb(loop_type, search_space_pdb_id, results, save_path)
            search_space_pdb_id = search_space_id.split("_")[1]
            results = {} # results: {(query_id,search_space_id) : [dq,discrepancy]}
            if load_saved_searches:
                timer_data = myTimer("Loading files")
                results = load_previous_search_results_one_pdb(loop_type, search_space_pdb_id, save_path)
            new_results = False

        for query_id in sorted(queries.keys()):
            if(query_id != search_space_id):

                # Check if this search already done
                if (query_id, search_space_id) not in results:
                    new_results = True
                    query_length = len(queries[query_id]['Q']['fullUnits'])
                    # query_length = queries[query_id]['Q']['numpositions']
                    query_bulge_length = len(queries[query_id]['loop_info']['bulged'])
                    is_single_bulged_candidate = query_bulge_length ==1 and query_length==5

                    search_space_length = len(search_spaces[search_space_id]['ifedata']['index_to_id'])
                    search_space_bulged_length = len(search_spaces[search_space_id]['loop_info']['bulged'])
                    is_single_bulged_search_space = search_space_bulged_length ==1 and search_space_length==5

                    if query_id[0:2] == 'IL' and is_single_bulged_candidate:
                        # treat single base bulge IL differently, to keep them together according to bulged base
                        if is_single_bulged_search_space:
                            Q, candidates, elapsed_time = FR3D_search(Q = queries[query_id]['Q'],
                                        ifedata = search_spaces[search_space_id]['ifedata'], ifename = search_space_id,
                                        timerData = timer_data)
                            if candidates:
                                temp_result = analyze_single_bulged_base_loops(candidates=candidates,query=queries[query_id],search_space=search_spaces[search_space_id])
                                if temp_result:
                                    temp_result[0]['match_type'] = "geometric"
                                    # temp_result[0]['match_type'] = "single base bulge"
                                    results[(query_id, search_space_id)] = temp_result
                            else:
                                results[(query_id, search_space_id)] = MISMATCHED_BULGE
                        else:
                            results[(query_id, search_space_id)] = MISMATCHED_BULGE
                    elif search_space_id[0:2] =="IL" and is_single_bulged_search_space:
                        results[(query_id, search_space_id)] = MISMATCHED_BULGE

                    else:
                        if(query_length - query_bulge_length <= search_space_length):
                            # Check flanking bp for IL and larger to make sure the match is plausible
                            if loop_type != 'HL':
                                Q,candidates, elapsed_time = FR3D_search(Q=flanking_bp_queries[query_id]['Q'],ifedata = search_spaces[search_space_id]['ifedata'], ifename = search_space_id,
                                timerData = timer_data)
                            else:
                                candidates = []

                            if loop_type == "HL" or candidates:
                                # do a full search for query inside of search_space
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
                                            # if temp_result:
                                            dq = temp_result[0]['dq']
                                            if len(dq) == 0:
                                                temp_result[0]['match_type'] = "geometric"
                                                results[(query_id,search_space_id)] = temp_result
                                            else:
                                                results[(query_id,search_space_id)] = dq[0]  # save one disqualification code
                                            # else:
                                            #     results[(query_id, search_space_id)] = UNMATCHED_BASEPAIR
                                        else:
                                            results[(query_id, search_space_id)] = UNMATCHED_BASEPAIR
                                    else:
                                        results[(query_id, search_space_id)] = CONFLICTING_BASEPAIRS_AND_STACKS

                                    if len(temp_result) > 0:
                                        print('Found %s inside %s in %5d ways, and %5d remain after filtering, the size of results is %s' % (query_id,search_space_id,len(candidates),len(temp_result), len(results)))

                                else: # no candidates from FR3D search
                                    #results[(query_id, search_space)] = [{'dq': 0, 'discrepancy': 99}] # hoping disqualification code of 0 makes sense for NO MATCH
                                    #Put 0 instead of list
                                    results[(query_id, search_space_id)] = NO_CANDIDATES
                                # after candidate sorting, we will no longer use Find_lowest()
                                # run test_motif_atlas_code::check_interaction() after EACH FR3D search Adam
                            else:
                                results[(query_id, search_space_id)] =  FLANKING_BP_CONFLICT
                        else:
                            results[(query_id, search_space_id)] =  SEARCH_SPACE_CONFLICT

                if isinstance(results[(query_id, search_space_id)],int):
                    search_results[(query_id,search_space_id)] = [{'dq': [NO_CANDIDATES], 'discrepancy': 99}]
                else:
                    search_results[(query_id,search_space_id)] = results[(query_id,search_space_id)]

        # moved up top. if file save issues, switch back to this block
        # if new_results:
        #     timer_data = myTimer('Saving search results according to search space pdb id')
        #     save_search_results_one_pdb(loop_type, search_space_pdb_id, results, save_path)
        #     new_results = False

    # save the results for the last PDB, which may be the only PDB
    if len(results) > 0:
        print('Saving %d results for %s and %s in %s' % (len(results), loop_type, search_space_pdb_id, save_path))
        save_search_results_one_pdb(loop_type, search_space_pdb_id, results, save_path)

    return(search_results)


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


def cluster_loops_hierarchical(disc_m,cluster_method):
    """
    Use hierarchical clustering to find tight groups of loops
    """

    N = len(disc_m)

    # get a list of mergers between individuals and clusters
    Z = linkage(squareform(disc_m), cluster_method.replace("hierarchical_",""))

    # track the different groups as they are merged together
    group_num_to_list = {}   # the members of the group, but emptied once merged

    for i in range(0,N):
        # the first groups are singletons
        group_num_to_list[i] = [i]

    # loop over the mergers and recognize additional groups
    for merger in Z:
        if merger[2] < 10:
            print('Considering a merger at distance %8.4f' % merger[2])
            # single bulged bases can have discrepancy above 1
            # just avoid whatever large number indicates a disqualification
            a = int(merger[0])
            b = int(merger[1])
            i += 1
            group_num_to_list[i] = group_num_to_list[a] + group_num_to_list[b]
            # now that these are merged, empty them
            group_num_to_list[a] = []
            group_num_to_list[b] = []

    groups = [g for g in group_num_to_list.values() if len(g) > 0]

    groups = sorted(groups, key = lambda x : len(x), reverse = True)

    for i,g in enumerate(groups):
        if len(g) > 1:
            print('Group %3d has size %3d' % (i,len(g)))

    print('Found %4d groups' % len(groups))

    return groups


def analyze_single_bulged_base_loops(candidates,query, search_space):
    not_rejected = []
    for candidate in candidates:
        DQ_code = []
        query_bulge,search_space_bulge = get_bulge(query,search_space)
        if query_bulge != search_space_bulge:
            DQ_code.append(MISMATCHED_BULGE)

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
    return(temp_result)


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






def search_for_centriod_loop(loop_ids,loop_ids_to_discrepancy, search_results,queries):
    min_total_disc = maxsize
    centroid_id = ""
    dict_for_units_sets = {}
    for oneLoop in loop_ids:
        oneLoopUnits = set(queries[oneLoop]["Q"]["unitID"]) # write all unit ids in oneLoop
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
        # if (total_disc < min_total_disc) and (len(oneLoopUnits) >= numberOfAlignment):
        if len(oneLoopUnits) >= largestNumberOfAlignment and (total_disc < min_total_disc):
            min_total_disc = total_disc
            largestNumberOfAlignment = len(oneLoopUnits)
            centroid_id = oneLoop
            dict_for_units_sets[oneLoop] = []
            for unit in queries[oneLoop]["Q"]["unitID"]:
                if unit in oneLoopUnits:
                    dict_for_units_sets[oneLoop].append(unit)
    return centroid_id, dict_for_units_sets

def align_candidate_units_to_centroid_units(real_ordered_unit_ids, centroid_unit_ids, candidate_unit_ids, centroid_to_candidate, core_positions_candidate_id):
    for centroid_unit_id in real_ordered_unit_ids:
        centroid_unit_id_index = centroid_unit_ids.index(centroid_unit_id)
        candidate_unit_id = candidate_unit_ids[centroid_unit_id_index]
        centroid_to_candidate[centroid_unit_id] = candidate_unit_id
        core_positions_candidate_id.append(centroid_to_candidate[centroid_unit_id])


def align_nts(loop_ids, loop_ids_to_discrepancy, search_results, queries):
    """
    Identify core positions in the given cluster and align all instances to that
    loop_ids is a list of text strings
    loop_ids_to_discrepancy is a
    """
    core_positions = defaultdict()

    # decide centroid loop id:
    centroid_id, dict_for_units_sets = search_for_centriod_loop(loop_ids, loop_ids_to_discrepancy, search_results, queries)

    # put the longer strand first
    dict_for_seperate_strand = {}

    for strand_index in range(len(queries[centroid_id]["loop_info"]["strand"])):
        dict_for_seperate_strand[strand_index] = []
    #need OneLoopUnits around here

    # add a condition for skipping HLs ? really need it ?
    # for strand in queries[centroid_id]["loop_info"]["strand"]:
    #     strand_name = queries[centroid_id]["loop_info"]["strand"].index(strand)
    #     for unit_id in strand:
    #         if unit_id in dict_for_units_sets[centroid_id]:
    #             dict_for_seperate_strand[strand_name].append(unit_id)

    for unit_id in dict_for_units_sets[centroid_id]:
        for strand_index, strand in enumerate(queries[centroid_id]["loop_info"]["strand"]):
            if unit_id in strand:
                dict_for_seperate_strand[strand_index].append(unit_id)

    # reorder nucleotides?
    max_key = max(dict_for_seperate_strand, key=lambda k: len(dict_for_seperate_strand[k]))
    real_ordered_unit_ids = []
    if max_key != 0:
        original_order = list(dict_for_seperate_strand.keys())
        new_order = original_order[max_key:] + original_order[:max_key]
        sorted_dict_for_seperate_strand = OrderedDict((k, dict_for_seperate_strand[k]) for k in new_order)
        for strand_name, unit_ids in sorted_dict_for_seperate_strand.items():
            real_ordered_unit_ids.extend(unit_ids)
    else:
        # original_order = list(dict_for_seperate_strand.keys())
        # new_order = original_order[max_key:] + original_order[:max_key]
        for strand_name, unit_ids in dict_for_seperate_strand.items():
            real_ordered_unit_ids.extend(unit_ids)
    core_positions[centroid_id] = real_ordered_unit_ids

    # align candidate loops with centroid loop
    for candidate_id in loop_ids:
        if candidate_id != centroid_id:
            centroid_to_candidate = {} # for target loop
            sr = search_results[(centroid_id,candidate_id)][0]
            if "target_unit_ids" in sr:
                candidate_unit_ids = sr["target_unit_ids"]
                centroid_unit_ids = sr["query_unit_ids"]
            else:
                sr = search_results[(candidate_id,centroid_id)][0]
                candidate_unit_ids = sr["query_unit_ids"]
                centroid_unit_ids = sr["target_unit_ids"]
            core_positions[candidate_id] = []
            align_candidate_units_to_centroid_units(real_ordered_unit_ids, centroid_unit_ids, candidate_unit_ids, centroid_to_candidate, core_positions[candidate_id])

    return core_positions, centroid_id


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


def save_motif_atlas_already_setup(release, atlas_queries, atlas_search_spaces, atlas_flanking_bp_queries):
    '''
    this saves some complicated dictionaries to avoid constant reprocessing
    '''
    if not os.path.exists(DATAPATHATLAS):
        os.mkdir(DATAPATHATLAS)

    pickle.dump(obj = atlas_queries, file = open(os.path.join(DATAPATHATLAS,
        "%s_atlas_queries.pickle" % release), "wb"), protocol = 2)
    pickle.dump(obj = atlas_search_spaces, file = open(os.path.join(DATAPATHATLAS,
        "%s_atlas_search_spaces.pickle" % release), "wb"), protocol = 2)
    pickle.dump(obj = atlas_flanking_bp_queries, file = open(os.path.join(DATAPATHATLAS,
        "%s_atlas_flanking_bp_queries.pickle" % release), "wb"), protocol = 2)
    return


def load_motif_atlas_already_setup(release):
    '''
    this is the counterpart to the above function "save_motif_atlas_already_setup"
    '''
    atlas_queries = pickle.load(open(os.path.join(DATAPATHATLAS,
        "%s_atlas_queries.pickle" % release), "rb"))
    atlas_search_spaces = pickle.load(open(os.path.join(DATAPATHATLAS,
        "%s_atlas_search_spaces.pickle" % release), "rb"))
    atlas_flanking_bp_queries = pickle.load(open(os.path.join(DATAPATHATLAS,
        "%s_atlas_flanking_bp_queries.pickle" % release), "rb"))

    return atlas_queries, atlas_search_spaces, atlas_flanking_bp_queries


def main(loops = None, pair_to_interaction_list = None, loops_and_strands = None, output_dir = './'):
    # loops_and_strands=load_loops_and_strands()

    queries = {}
    search_spaces ={}
    flanking_bp_queries = {}
    flanking_bp_queries2 = {}

    load_all_previous_results = True
    load_all_previous_results = False

    # to search faster, run two instances of this program with True/False on next line
    reversed_search = True
    reversed_search = False

    cluster_method = 'clique_maximum'
    cluster_method = 'clique_average'
    cluster_method = 'hierarchical_complete'
    cluster_method = 'hierarchical_average'

    ratio = 0.9          # the ratio that determines what range of clique sizes to compare

    timer_data = myTimer("start")
    timer_data = myTimer('Set up loops')

    # focus on specific cases
    loops_of_interest = [loop_id for loop_id,alignment in loops_and_strands.items() if len(alignment) == 5]
    loops_of_interest = ['IL_7RQB_111','IL_5J7L_302','IL_4WF9_002', 'IL_4V88_409', 'IL_4V9F_002']
    loops_of_interest = ["IL_7A0S_074","IL_7RQB_079"]   # they should match well, but do they?
    loops_of_interest = ["J3_4V88_034","J3_4LFB_005", "J3_5J7L_007"]
    loops_of_interest = ["J3_7RQB_004","J3_7A0S_004", "J3_4V9F_002", "J3_5J7L_038", "J3_5TBW_004", "J3_4WF9_002"]
    loops_of_interest = ["IL_5J7L_297","IL_5J7L_303"]
    loops_of_interest = ["IL_6UFM_002","IL_6UFG_002","IL_6UFH_002"]
    load_saved_searches = False    # force the program to re-compute the searches, good for testing

    # analyze all loops loaded
    loops_of_interest = [loop_id for loop_id,alignment in loops_and_strands.items()]
    load_saved_searches = False    # force the program to re-compute the searches, good for testing
    load_saved_searches = True     # load results of previous searches if available, faster

    # Xinyu's testing
    # if 'xchu' in os.getcwd():
    #     Chu_test = True
    # else:
    #     Chu_test = False
    # Chu_test = False
    # if Chu_test:
    #     loops_of_interest = ["IL_1MFQ_001","IL_1L9A_001","IL_1KUQ_003", "IL_6CZR_141", 'IL_4LFB_029',"IL_1G1X_006", 'IL_1G1X_003', 'IL_5J7L_033', 'IL_4TS2_005', 'IL_4KZD_004', 'IL_5TBW_043']
    #     # loops_of_interest = [loop_id for loop_id,alignment in loops_and_strands.items()]
    #     loops_of_interest=['IL_5TBW_077','IL_4V9F_062','IL_4WF9_065','IL_5J7L_307','IL_7A0S_062','IL_7RQB_067','IL_7A0S_044','IL_4WF9_047','IL_5VCI_005','IL_2O3X_001','IL_2ET8_001','IL_2O3X_003','IL_7A0S_092',]
    #     load_saved_searches = True

    # if there is user input, it's more important than the above tests
    if loops:
        loops_of_interest = loops

    loop_type = set_loop_type(loops_and_strands)

    release = '3.57'

    # load motif atlas release
    motif_atlas_release = load_motif_atlas_release(loop_type,release)
    motif_id_to_loops,loop_id_to_motif_id = motif_ids(motif_atlas_release,loops_of_interest)

    # load loop annotations
    loop_id_to_annotation = defaultdict(str)
    annotation_filename = os.path.join(DATAPATHLOOPS, "loop_annotations.txt")
    with open(annotation_filename,"r") as file:
        for line in file:
            data = line.split("\t")
            loop_id_to_annotation[data[0]] = data[1]

    # evaluate how well the Matlab groups do, with loop annotations
    matlab_groups = []
    for motif_id, loops in motif_id_to_loops.items():
        matlab_groups.append(loops)
    print('Evaluation of Matlab grouping:')
    evaluate_grouping(matlab_groups,loop_id_to_annotation)


    if load_all_previous_results:
        # this may no longer be needed since the rest of the code is much faster now - CLZ 2023-06-22
        search_results,all_loop_ids = load_all_search_results(loop_type)
        queries = load_all_queries(loop_type)

        #reordered_motif_groups = load_all_clusters(loop_type)
        #reordered_dist_matrices = load_dist_matrices(loop_type)
        #reordered_motif_groups,reordered_dist_matrices = order_groups_by_similarity(reordered_motif_groups,search_results)
        #validate_clusters(search_results,loops_of_interest,reordered_motif_groups)
        show_discrepancy(search_results, ['IL_6CZR_120'],['IL_1QBP_002','IL_1IHA_001','IL_1ICG_001'])

    else:
        print("Testing %d %s loops of interest" % (len(loops_of_interest),loop_type))

        # process a raw list of loops and their nucleotides
        # load FR3D annotations of pairwise interactions
        loops, pair_to_interaction_list = startup_list_of_dictionaries(loops_and_strands,loops_of_interest)

        # sort the list of loops by loop id
        loops = sorted(loops, key = lambda x : x['loop_id'])

        timer_data = myTimer("Set up queries and ifedata",timer_data)
        print("Set up queries and search space data")

        # for each loop, set it up as a query Q and as a search space ifedata
        all_queries, all_search_spaces, flanking_bp_queries = create_queries_and_search_spaces(loops,loops_of_interest)

        queries = {}
        search_spaces = {}
        for loop_id in loops_of_interest:
            if loop_id in all_queries:
                queries[loop_id] = all_queries[loop_id]
                search_spaces[loop_id] = all_search_spaces[loop_id]

        save_all_queries(loop_type,queries)

        print("Start all against all searches")
        timer_data = myTimer("All against all searches",timer_data)
        motif_result_path = MOTIFPATH
        search_results = all_against_all_searches(loop_type,queries,search_spaces,flanking_bp_queries,load_saved_searches,reversed_search, motif_result_path)

    # index the loop ids in a list, so you can refer to them by number
    all_loop_ids = list(queries.keys())

    print("Found %d loops with full coordinate data" % len(all_loop_ids))


    timer_data = myTimer("Calculate distance matrices",timer_data)
    print("Calculating distance matrices")

    disc_m, MM, dq_matrix = calculate_distance_matrices(search_results, all_loop_ids)


    # put literally all of the loops into one ordering; can be a useful diagnostic
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


    # map loop id to index
    loop_id_to_index = {}
    for i, loop_id in enumerate(all_loop_ids):
        loop_id_to_index[loop_id] = i


    if len(loops_of_interest) < 10:
        print("Too few loops of interest, not doing additional diagnostics")
        return


    # diagnostic over Matlab motif groups
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
        motif_groups = cluster_loops_hierarchical(disc_m,cluster_method)

    # find maximum distance within each group
    group_max_distance = {}
    for i, motif_group in enumerate(motif_groups):
        maxd = 0
        for a in motif_group:
            for b in motif_group:
                maxd = max(maxd,disc_m[a][b])
        group_max_distance[i] = maxd


    # find minimum distance between groups
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

    # list motif group by loop ids
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

    # evaluate how well the new clusters do, with loop annotations
    print("Evaluation of how well the new clusters work")
    evaluate_grouping(motif_groups_by_loop_id,loop_id_to_annotation)

    print("Mapped %4d loop ids to motif groups, out of %4d loop ids" % (len(loop_id_to_group_number.keys()),len(all_loop_ids)))
    missed_loop_ids = set(all_loop_ids) - set(loop_id_to_group_number.keys())
    if len(missed_loop_ids) > 0:
        print('There are %3d missing loop ids:' % len(missed_loop_ids))
        print(sorted(missed_loop_ids))

    # identify groups a loop could be in (including its own) and get max and mean discrepancy
    loop_id_to_matching_groups = defaultdict(dict)
    for i in range(0,len(all_loop_ids)):
        loop_id = all_loop_ids[i]
        group_num_to_discrepancy = defaultdict(list)

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

    # re-order each group, find core positions, write HTML file for inspection
    for num, motif_group in enumerate(motif_groups):

        # re-order distance matrices and change indexing from numeric to loop ids
        reordered_motif_group, reordered_dist_matrix = order_group_by_similarity(motif_group,disc_m)

        # index motif group and discrepancy matrix by loop ids
        reordered_loop_ids = [all_loop_ids[x] for x in reordered_motif_group]

        loop_ids_to_discrepancy = {}
        for i, loop_i_id in enumerate(reordered_loop_ids):
            loop_ids_to_discrepancy[loop_i_id] = {}
            for j, loop_j_id  in enumerate(reordered_loop_ids):
                loop_ids_to_discrepancy[loop_i_id][loop_j_id] = reordered_dist_matrix[i][j]

        # save .pickle files in DATAPATH/new_results dictory
        # Slightly faster to save all search results and then re-load, but savefile is > 100 MB
        #save_all_search_results(loop_type,search_results,all_loop_ids)
        #save_dist_matrices(loop_type,reordered_dist_matrices)
        #save_all_clusters(loop_type,reordered_motif_groups,cluster_method)
        #show_discrepancy(search_results, ['IL_3IGI_011'],['IL_4V9F_007'])

        # identify core positions in each group
        # Why?  The instances in a cluster might not all have the same number of nucleotides
        # It's common for a loop in one structure, one organism to have one or more "bulged out" nucleotides
        # The FR3D search technique we use is robust to those bulged out nucleotides
        # The FR3D search basically gives pairwise alignments between loop instances
        # But for a motif group, we need a complete, single alignment, not a bunch of pairwise alignments
        # We need to identify "core" and "non-core" nucleotides

        #save_all_search_results(search_results,all_loop_ids) #MAYBE LATER
                
        # Xinyu, at this point you could write out some .pickle files with the data you need
        # and then read those files and work on consensus basepairs and whatever
        core_positions, centroid = align_nts(reordered_loop_ids, loop_ids_to_discrepancy, search_results, queries)
        
        matrixForConsensusInteractions = get_consensus_interactions(core_positions)
        # writeHTMLOutput(num,loop_type,cluster_method,search_results,motif_groups,core_positions,reordered_loop_ids,queries,loop_ids_to_discrepancy,loop_id_to_matching_groups,loop_id_to_annotation,loop_id_to_group_number,group_max_distance,release)
        folder_name = output_dir + '/2ds'
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

        num = num+1
        varna_command_info = motif_to_varna(core_positions, matrixForConsensusInteractions, cluster_method)



        varna_command = 'java -cp bin/VARNAv3-7.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN "%s" -structureDBN "%s" -baseNum "#334455" -periodNum 1 -auxBPs "%s" -o %s/2ds/Group_%03d.svg' % (varna_command_info['sequenceDBN'], varna_command_info['structureDBN'], varna_command_info['auxBPs'], output_dir, num)
        os.system(varna_command)
        varna_command = 'java -cp bin/VARNAv3-7.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN "%s" -structureDBN "%s" -baseNum "#334455" -periodNum 1 -auxBPs "%s" -o %s/2ds/Group_%03d.png' % (varna_command_info['sequenceDBN'], varna_command_info['structureDBN'], varna_command_info['auxBPs'], output_dir, num)
        os.system(varna_command)
        with open(output_dir + '/' + 'varna_commands.txt', 'a+') as file:
            file.write(varna_command + '\n')

        # debugging problemtic varna commands:
        varna_image_path = '%s/2ds/Group_%03d.svg' % (output_dir, num)
        if not os.path.exists(varna_image_path):
            with open(output_dir + '/' + 'problemtic_varna_commands.txt', 'a+') as file:
                file.write(varna_command + '\n')
            first_varna_image_path = '/usr/local/pipeline/hub-core/pymotifs/motif_atlas/unkown_varna.png'
            shutil.copy(first_varna_image_path, varna_image_path)
        # bpSignature = find_bp_signature(queries, core_positions, centroid, matrixForConsensusInteractions)
        print('debug info: the centroid is %s, the actual units are: %s, the core_positions is %s, the matrix for consensus interactions are: %s' % (centroid, set(queries[centroid]["Q"]["unitID"]), core_positions, matrixForConsensusInteractions))
        writeCSVOutput(group_id = num, loop_ids = reordered_loop_ids, loop_ids_to_discrepancy = loop_ids_to_discrepancy, core_positions= core_positions, centroid_loop = centroid, matrixForConsensusInteractions = matrixForConsensusInteractions, queries=queries, output_dir = output_dir)
        # writeHTMLOutput(num,loops_of_interest, cluster_method,search_results,motif_groups,core_positions,reordered_loop_ids,queries,loop_ids_to_discrepancy,loop_id_to_matching_groups,loop_id_to_annotation,loop_id_to_group_number,group_max_distance,release, matrixForConsensusInteractions,varna_command_info,bpSignature)
        '''
        test
        '''


    # get statistics about agreement between Python clusters and Matlab motif groups
    # output is written to DATAPATH/new_results directory

    print("Remember to fix up the code in validate_clusters")
    # print("Note that align_nts will sometimes give different numbers of core nucleotides!")
    # print("Note that align_nts will sometimes completely leave out nucleotides from a loop, neither marked as core or bulged!")

    """
    timer_data = myTimer("Validate clusters",timer_data)
    validate_clusters(cluster_method,loop_type,search_results,loops_of_interest,reordered_motif_groups)
    """

    print(myTimer("summary"))


if __name__ == '__main__':
    from HL_10_24_loops_and_strands import loops_and_strands
    from IL_10_05_loops_and_strands import loops_and_strands
    from IL_3_57_loops_and_strands import loops_and_strands


    res = main(loops_and_strands = loops_and_strands)
