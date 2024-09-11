#from loops_and_strands2 import all_structures, loops, pair_to_interaction_list
import fr3d_interactions
from collections import defaultdict

all_bptypes = {'cWW', 'tWW', 'cWH', 'tWH', 'cWS', 'tWS', 'cHH', 'tHH', 'cHS', 'tHS',\
                   'cSS', 'tSS','cHW','tHW','cSW','tSW','cSH','tSH'}
near_bptypes = {'ntHS', 'ntSH', 'ntWH', 'ncHS', 'ncSH', 'ncWW', 'ncHH', 'ntWW', 'ntSS',\
                    'ncSS', 'ncWS', 'ntWS', 'ntHH', 'ntHW', 'ncSW', 'ncHW', 'ncWH', 'ntSW'}
stacks = {'s33','s35','s55','s53'}
near_stacks = {'ns35','ns55','ns33','ns53'}


def add_bulged_nucleotides(loops_lst, inter_lst):
    """
    Returns the list of dictionaries of loops that will now contain a list of
    the bulged nucleotide unit ids.
    {id: HL_4V9F_003, strands: [[4V9F|1|0|G|421,4V9F|1|0|G|422,4V9F|1|0|C|423]],bulged: [4V9F|1|0|G|422]}
    """
    non_bulged_interactions = all_bptypes | near_bptypes | stacks | near_stacks

    for loop in loops_lst:

        bulged_nts = []
        all_nts = [unit_id for sublist in loop['strand'] for unit_id in sublist]

        for unit1 in all_nts:
            unit1_interactions = set()

            for unit2 in all_nts:
                unit1_interactions.update(inter_lst[(unit1, unit2)])

            if not unit1_interactions.intersection(non_bulged_interactions):
                bulged_nts.append(unit1)
        loop['bulged'] = bulged_nts

    return loops_lst

# insert code to use loops list to make aAa searches!!! ----------------------------------

def annotate_pair_stack_conflicts(aAa_searches,inter_lst):
    """
    This method ensures that all of the aAa searches have another key, called dismilarity,
    where the value is either a string indicating a basepair or basestack conflict or the
    discrepancy if there is no such conflict.

    aAa_searches[('IL_4V9F_001','IL_4V9F_025')] =  {matched1: [matched unit ids in
    1st loop], matched2: [matched unit ids in 2nd loop],'disc':discrepancy, 'extra_nts1':
    [lst of extra nts in 1st loop],'extra_nts2':[lst of extra nts in 2nd loop]}

    """
    aAa_searches_tuples = aAa_searches.keys()
    for loop_tuple in aAa_searches_tuples:
        # if check_interactions returns something, which means there is a conflict
        temp = check_interactions(aAa_searches[loop_tuple],inter_lst)
        if temp:
            aAa_searches[loop_tuple]['dissimilarity'] = temp
        # if not, means there's no conflict, so disimilarity will just be the discrepancy value
        else:
            aAa_searches[loop_tuple]['dissimilarity'] = aAa_searches[loop_tuple]['disc']
    return aAa_searches


def check_interactions(pairwise_search,inter_lst):
    """
    Goes through the possible interactions made between loop1 and loop2 using 2 for loops
    and checks if there are conflicting base pairs or stacking.

    aAa_searches[('IL_4V9F_001','IL_4V9F_025')] =  {matched1: [matched unit ids in
    1st loop], matched2: [matched unit ids in 2nd loop],'disc':discrepancy, 'extra_nts1':
    [lst of extra nts in 1st loop],'extra_nts2':[lst of extra nts in 2nd loop]}
    """
    matched_units1 = pairwise_search['matched1'] # query
    matched_units2 = pairwise_search['matched2'] # target
    num_units = len(matched_units1)

    for i in range(num_units):
        for j in range(num_units):
            if i != j:
                inter_i = inter_lst[(matched_units1[i],matched_units1[j])]
                inter_j = inter_lst[(matched_units2[i],matched_units2[j])]

                if inter_i and inter_j and set(inter_i) != set(inter_j):

                    # if both interactions in the same aligned position of loop i and j are basepairs,
                    # but different/ conflicting.
                    if inter_i in all_bptypes and inter_j in all_bptypes and inter_i != inter_j:
                        return("basepair conflict: {}  {}".format(inter_i, inter_j))

                    # if interaction in loop i is a stack and interaction in loop j is a basepair,
                    # vice versa
                    if (inter_i in all_bptypes and inter_j in stacks) or (inter_i in stacks and inter_j in all_bptypes):
                        return("basestack conflict: {}  {}".format(inter_i, inter_j))
    return None


def check_aAa_extra_nts(aAa_searches,inter_lst):
    """
    Parameters: aAa_searches: dictionary of dictionaries, key is tuple of 2
    loops. aAa_searches[('IL_4V9F_001','IL_4V9F_025')] =  {matched1: [matched unit ids in
    1st loop], matched2: [matched unit ids in 2nd loop],'disc':discrepancy, 'extra_nts':
    [lst of extra nts],'core_nts':[lst of core nts]}

    """
    # write a method to look through all aAa searches and analyze extra nts using method below,
    # and change ['disimilarity'] value
    # add extra nts1 and extra nts2, delete core_nts
    # matched1 will be a list of lists
    aAa_searches_tuples = aAa_searches.keys()
    for loop_tuple in aAa_searches_tuples:
        #loop tuple is loop type, first 2 chars of loop ID
        temp = analyzeExtraNts(inter_lst,aAa_searches[loop_tuple]['extra_nts'],\
                               aAa_searches[loop_tuple]['core_nts'],loop_tuple[0][0:2])
        if temp:
            aAa_searches[loop_tuple]['dissimilarity'] = temp
        else:
            aAa_searches[loop_tuple]['dissimilarity'] = aAa_searches[loop_tuple]['disc']

    return aAa_searches

def analyzeExtraNts(inter_lst,extraNts_ids,coreNts,loop_type):
    """
    Parameters: All interactions list: inter_lst[(nt,nt)] gives list of
    interactions made between nt1 and nt2, extraNts_ids: list of the
    NT IDs of the extra nucleotides, coreNts: list of NT IDs of the core
    nucleotides
    This method checks if the extra nucleotides make bad interactions, and
    returns a penalty assigned if it does.
    """

    interactionsWithCore = [] #list of interactions extra nts make with core
    interactionsInExtra = [] #list of interactions extra nts make with other extra nts

    for extra_nt in extraNts_ids:
        for core_nt in coreNts:
            if (extra_nt,core_nt) in inter_lst.keys():
                interactionsWithCore.extend(inter_lst[(extra_nt,core_nt)])
        for extra_nt2 in extraNts_ids:
            if (extra_nt,extra_nt2) in inter_lst.keys():
                interactionsInExtra.extend(inter_lst[(extra_nt,extra_nt2)])

    interactionsWithCore = set(interactionsWithCore)
    interactionsInExtra = set(interactionsInExtra)
    stackCount = len(interactionsWithCore.intersection(stacks))
    bpInteractionsWithCore = interactionsWithCore.intersection(all_bptypes)
    bpInteractionsInExtra = interactionsInExtra.intersection(all_bptypes)

    # extra nucleotides make bps with any nucleotide in the motif
    if bpInteractionsWithCore:
        return("extra nucleotide makes {} interaction".format(bpInteractionsWithCore))
    elif bpInteractionsInExtra: # extra nts interacting on eachother
        return("extra nucleotide makes {} interaction".format(bpInteractionsInExtra))

    # extra nucleotides make near bps with the core nucleotides
    elif interactionsWithCore.intersection(near_bptypes):
        return("extra nucleotide makes near basepair with core nucleotide")

    # extra nucleotides make more than 1 stacking interaction with coreNts
    elif stackCount > 1:
        return("extra nucleotide makes {} stacking interactions with coreNts".format(stackCount))

    # extra nucleotides in hairpins stack on coreNts
    elif loop_type == 'HL' and interactionsWithCore.intersection(stacks):
        return("extra nucleotide in hairpin stacks on core nucleotide")

    # extra nucleotides in hairpins stack on each other
    elif loop_type == 'HL' and interactionsInExtra.intersection(stacks):
        return("extra nucleotide in hairpin stacks on each other")


def pretend_aAa_searches(loops,pair_to_interaction_list):
    """
    Makes test aAa searches using loops that have both strands with length 5, and making the
    default discrepancy 0.2
    aAa_searches[('IL_4V9F_001','IL_4V9F_025')] =  {
                                       'matched1': list of lists of matched unit ids in 1st loop,
                                       'matched2': list of lists of matched unit ids in 2nd loop,
                                       'disc': discrepancy,
                                       'extra_nts': [lst of extra nts],
                                       'core_nts':[lst of core nts]
    }
    """
    aAa_searches = {}
    loop_ids = []

    loops = add_bulged_nucleotides(loops, pair_to_interaction_list)

    for i in range(len(loops)):
        for j in range(len(loops)):
            # length of strands
            strand1_len_i = len(loops[i]['strand'][0])
            strand1_len_j = len(loops[j]['strand'][0])
            strand2_len_i = len(loops[i]['strand'][1])
            strand2_len_j = len(loops[j]['strand'][1])

            if strand1_len_i == 5 and strand1_len_j == 5 and strand2_len_i== 5 and strand2_len_j == 5:
                aAa_searches[(loops[i]['loop_id'], loops[j]['loop_id'])] = {
                    'matched1':loops[i]['strand'][0]+loops[i]['strand'][1],
                    'matched2':loops[j]['strand'][0]+loops[j]['strand'][1],
                    'disc': 0.2,
                    'extra_nts1':loops[i]['bulged'],
                    'extra_nts2':loops[j]['bulged']
                }
                if loops[i]['loop_id'] not in loop_ids:
                    loop_ids.append(loops[i]['loop_id'])
                if loops[j]['loop_id'] not in loop_ids:
                    loop_ids.append(loops[j]['loop_id'])

    return (aAa_searches,loop_ids)

def main():
    pair_stack_list = defaultdict(list)

    # pair_stack_list is a dictionary of pair and stacking interactions
    for pdb_id in all_structures:
        pair_stack_list.update(fr3d_interactions.get_fr3d_pair_to_basepair_or_stack(pdb_id))

    (aAa_searches,loop_ids) = pretend_aAa_searches(loops,pair_to_interaction_list)
    aAa_searches = annotate_pair_stack_conflicts(aAa_searches,pair_stack_list)
    return()
