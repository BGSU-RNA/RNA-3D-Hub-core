import os
import pickle
import sys
import urllib
from collections import defaultdict
sys.path.append('search')
from fr3d_configuration import DATAPATHPAIRS

def get_fr3d_interaction_to_pair_list(pdb_id):
    """
    This method downloads FR3D annotations if necessary,
    then reads the .pickle data file and returns a dictionary
    which has interaction types like cWW and 5BPh as keys
    and whose values are lists of *triples* of unit ids which
    make the given interaction and the number of nested cWW
    basepairs which the interaction crosses.
    The third element of the triple is the crossing number / range.
    """

    pdb_id = pdb_id.upper()

    fr3d_interaction_file = "%s_RNA_pairs.pickle" % pdb_id
    fr3d_interaction_path_file = os.path.join(DATAPATHPAIRS, fr3d_interaction_file)

    if not os.path.exists(DATAPATHPAIRS):
        os.mkdir(DATAPATHPAIRS)

    if not os.path.isfile(fr3d_interaction_path_file):
        print("Downloading %s" % fr3d_interaction_file)
        url = "http://rna.bgsu.edu/pairs/%s" % fr3d_interaction_file
        if sys.version_info[0] < 3:
            urllib.urlretrieve(url, fr3d_interaction_path_file)  # python 2
        else:
            urllib.request.urlretrieve(url, fr3d_interaction_path_file)  # python 3

        if os.path.isfile(fr3d_interaction_path_file):
            try:
                with open(fr3d_interaction_path_file, 'rb') as opener1:
                    fr3d_interaction_to_pair_list = pickle.load(opener1)
            except:
                print("Unable to download %s" % fr3d_interaction_file)
                try:
                    os.remove(fr3d_interaction_path_file)
                except:
                    print("Unable to delete %s" % fr3d_interaction_file)
                return None

    with open(fr3d_interaction_path_file, 'rb') as opener1:
        fr3d_interaction_to_pair_list = pickle.load(opener1)

    return fr3d_interaction_to_pair_list

def get_fr3d_pair_to_basepair_type(pdb_id, near = True):
    """
    Returns a dictionary mapping pairs of unit ids to the
    basepair or near basepair they make, if any.
    Only the 12 Leontis-Westhof basepairs are included here.
    If no basepair is made, a blank string is returned.
    """

    all_bptypes= ['cWW', 'tWW', 'cWH', 'tWH', 'cWS', 'tWS', 'cHH', 'tHH', 'cHS', 'tHS', 'cSS', 'tSS']

    fr3d_interaction_to_pair_list = get_fr3d_interaction_to_pair_list(pdb_id)

    fr3d_pair_to_basepair_type = defaultdict(str)

    for bp_type in all_bptypes:
        # add the pair in the specified order, like cWS
        for (u1,u2,c) in fr3d_interaction_to_pair_list[bp_type]:
            fr3d_pair_to_basepair_type[(u1,u2)]=bp_type
        # add the pair in reversed order, like cSW
        reversed_bp_type=bp_type[0]+bp_type[2]+bp_type[1]
        for (u1,u2,c)  in fr3d_interaction_to_pair_list[bp_type]:
            fr3d_pair_to_basepair_type[(u2,u1)]=reversed_bp_type
        if near:
            # add near basepairs, like ncWS
            near_bp_type="n"+bp_type
            for (u1,u2,c) in fr3d_interaction_to_pair_list[near_bp_type]:
                fr3d_pair_to_basepair_type[(u1,u2)]=near_bp_type
            # add reversed near basepairs, like ncSW
            near_reversed_bp_type="n"+reversed_bp_type
            for(u1,u2,c) in fr3d_interaction_to_pair_list[near_bp_type]:
                fr3d_pair_to_basepair_type[(u2,u1)]=near_reversed_bp_type

    return fr3d_pair_to_basepair_type

def get_fr3d_pair_to_basepair_or_stack(pdb_id, near_basepair = True, near_stack = True):
    """
    Returns a dictionary mapping pairs of unit ids to the
    basepair or near basepair they make, if any.
    Only the 12 Leontis-Westhof basepairs are included here.
    If no basepair is made, a blank string is returned.
    """

    all_bptypes= ['cWW', 'tWW', 'cWH', 'tWH', 'cWS', 'tWS', 'cHH', 'tHH', 'cHS', 'tHS', 'cSS', 'tSS']
    all_stacks= ['s35','s55','s33']
    fr3d_interaction_to_pair_list = get_fr3d_interaction_to_pair_list(pdb_id)

    fr3d_pair_to_basepair_stack = defaultdict(str)

    for bp_type in all_bptypes:
        # add the pair in the specified order, like cWS
        for (u1,u2,c) in fr3d_interaction_to_pair_list[bp_type]:
            fr3d_pair_to_basepair_stack[(u1,u2)]=bp_type
        # add the pair in reversed order, like cSW
        reversed_bp_type=bp_type[0]+bp_type[2]+bp_type[1]
        for (u1,u2,c)  in fr3d_interaction_to_pair_list[bp_type]:
            fr3d_pair_to_basepair_stack[(u2,u1)]=reversed_bp_type
        if near_basepair:
            # add near basepairs, like ncWS
            near_bp_type="n"+bp_type
            for (u1,u2,c) in fr3d_interaction_to_pair_list[near_bp_type]:
                fr3d_pair_to_basepair_stack[(u1,u2)]=near_bp_type
            # add reversed near basepairs, like ncSW
            near_reversed_bp_type="n"+reversed_bp_type
            for(u1,u2,c) in fr3d_interaction_to_pair_list[near_bp_type]:
                fr3d_pair_to_basepair_stack[(u2,u1)]=near_reversed_bp_type
    
    for stack in all_stacks:
        # add the stack in the specified order
        for (u1,u2,c) in fr3d_interaction_to_pair_list[stack]:
            fr3d_pair_to_basepair_stack[(u1,u2)]=stack
        # add the pair in reversed order, like cSW
        reversed_stack=stack[0]+stack[2]+stack[1]
        for (u1,u2,c)  in fr3d_interaction_to_pair_list[stack]:
            fr3d_pair_to_basepair_stack[(u2,u1)]=reversed_stack
        if near_stack:
            # add near stacks
            near_stack="n"+stack
            for (u1,u2,c) in fr3d_interaction_to_pair_list[near_stack]:
                fr3d_pair_to_basepair_stack[(u1,u2)]=near_stack
            # add reversed near stacks
            near_reversed_stack="n"+reversed_stack
            for(u1,u2,c) in fr3d_interaction_to_pair_list[near_stack]:
                fr3d_pair_to_basepair_stack[(u2,u1)]=near_stack
    return fr3d_pair_to_basepair_stack

def get_fr3d_pair_to_interaction_list(pdb_id, near = True):
    """
    Returns a dictionary mapping pairs of unit ids to a list
    of all interactions they make, including base pairing,
    stacking, base-phosphate, and base-ribose, and including
    near versions of these interactions.
    Because a given pair can make both a basepair and a base-
    backbone interaction, the interactions are returned as a
    list.
    Note that 0BPh can be made as a self interaction between
    a nucleotide and itself.
    If the pair makes no interaction, an empty list is returned.

    """

    fr3d_interaction_to_pair_list = get_fr3d_interaction_to_pair_list(pdb_id)

    fr3d_pair_to_interaction_list = defaultdict(list)
    fr3d_pair_to_crossing_number = defaultdict(list)

    if not fr3d_interaction_to_pair_list: # new safety 09/21/2023
        return {}, {}

    for interaction_type in fr3d_interaction_to_pair_list.keys():
        for (u1,u2,c) in fr3d_interaction_to_pair_list[interaction_type]:
            fr3d_pair_to_interaction_list[(u1,u2)].append(interaction_type)
            fr3d_pair_to_crossing_number[(u1,u2)].append(c)


    return fr3d_pair_to_interaction_list, fr3d_pair_to_crossing_number

def demonstrate():
    # examples of using the methods above

    fr3d_interaction_to_pair_list = get_fr3d_interaction_to_pair_list("4V9F")
    print("Pairs making a cHH interaction and their crossing number:")
    print(fr3d_interaction_to_pair_list["cHH"])

    fr3d_pair_to_basepair_type = get_fr3d_pair_to_basepair_type("4V9F",False)
    print("Interaction between 4V9F|1|0|A|1742 and 4V9F|1|0|G|2033 is %s" % fr3d_pair_to_basepair_type[('4V9F|1|0|A|1742', '4V9F|1|0|G|2033')])
    print("Interaction between 1ABC|1|0|A|1742 and 1ABC|1|0|G|2033 is %s" % fr3d_pair_to_basepair_type[('1ABC|1|0|A|1742', '1ABC|1|0|G|2033')])

    fr3d_pair_to_interaction_list = get_fr3d_pair_to_interaction_list("4V9F")

    print("Pairs of nucleotides that make 3 or more interactions:")
    for a,b in fr3d_pair_to_interaction_list.items():
        if len(b) > 2:
            print(a,b)

    print("Interactions between 4V9F|1|0|G|64 and 4V9F|1|0|A|70 are %s" % fr3d_pair_to_interaction_list[('4V9F|1|0|G|64', '4V9F|1|0|A|70')])
    print("Interactions between 1ABC|1|0|G|64 and 1ABC|1|0|A|70 are %s" % fr3d_pair_to_interaction_list[('1ABC|1|0|G|64', '1ABC|1|0|A|70')])

#demonstrate()