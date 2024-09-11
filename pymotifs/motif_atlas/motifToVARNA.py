from __future__ import division
from asyncore import loop
from code import interact
from logging import logProcesses
from operator import index
from webbrowser import get

from collections import defaultdict

import fr3d_interactions
from fr3d_interactions import get_fr3d_interaction_to_pair_list, get_fr3d_pair_to_interaction_list
# import urllib.request
import sys
if sys.version_info[0] < 3:
    import urllib
else:
    import urllib.request
import numpy as np
# import pandas as pd
import os


def get_pdb_id(loop_id):
    pdb_id = loop_id.split('_')[1]
    return pdb_id



def get_consensus_interactions(core_positions):
    # pairTypes = ['glycosidicBondOrientation','chiDegree','pairsStacks','BPh','BR','sO','crossingNumber']
    pairTypes = ['basepair','pairsStacks','BPh','BR','others'] # change names stacks
    all_bptypes = ['ncHW', 'cSS', 'cWW', 'cWH', 'cHH', 'cHS', 'ntHW', 'ntHH', 'cSH', 'ntWS', 'ncSS', 'ncWH', 'ntWH', 'cHW', 'cWS', 'tHS', 'ntSH', 'tSH', 'tSW', 'ncHH', 'ncSH', 'tHW', 'ntSW', 'tHH', 'cSW', 'ncHS', 'ncWW', 'ntWW', 'tWS', 'ntSS', 'ntHS', 'tWW', 'ncSW', 'tSS', 'tWH', 'ncWS']
    all_stacks= ['s35','s55','s33','s53','ns35', 'ns55', 'ns33']

    all_unit_pairs ={}
    all_unit_pairs['basepair'] = set()
    all_unit_pairs['pairsStacks'] = set()
    all_unit_pairs['BR'] = set()
    all_unit_pairs['BPh'] = set()
    all_unit_pairs['others'] = set()

    matrixForConsensusInteractions = defaultdict()
    for loop_id, nts in core_positions.items():
        pdb_id = get_pdb_id(loop_id)
        fr3d_pair_to_interaction_list, fr3d_pair_to_crossing_number = get_fr3d_pair_to_interaction_list(pdb_id)
        matrixForConsensusInteractions[loop_id] = {pairType: {} for pairType in pairTypes}
        for oneUnit in range(len(nts)):
            for anotherUnit in range(oneUnit+1, len(nts)):
                pair = (nts[oneUnit], nts[anotherUnit])
                consensus_interaction = fr3d_pair_to_interaction_list[pair]
                if len(consensus_interaction) >= 1:
                    unitPairName = str(str(oneUnit + 1) + ' - ' + str(anotherUnit + 1))
                    for i in consensus_interaction:
                        if i in all_bptypes:
                            matrixForConsensusInteractions[loop_id]['basepair'][unitPairName] = i
                            all_unit_pairs['basepair'].add(unitPairName)
                        elif i in all_stacks:
                            matrixForConsensusInteractions[loop_id]['pairsStacks'][unitPairName] = i
                            all_unit_pairs['pairsStacks'].add(unitPairName)
                        elif 'BR' in i:
                            matrixForConsensusInteractions[loop_id]['BR'][unitPairName] = i
                            all_unit_pairs['BR'].add(unitPairName)
                        elif 'BPh' in i:
                            matrixForConsensusInteractions[loop_id]['BPh'][unitPairName] = i
                            all_unit_pairs['BPh'].add(unitPairName)
                        else:
                            matrixForConsensusInteractions[loop_id]['others'][unitPairName] = i
                            all_unit_pairs['others'].add(unitPairName)
                        
    # fill up all NULL value
    for pairType in pairTypes:
        for loop_id, nts in core_positions.items():
            for unitPairName in all_unit_pairs[pairType]:
                if unitPairName not in matrixForConsensusInteractions[loop_id][pairType]:
                    matrixForConsensusInteractions[loop_id][pairType][unitPairName] = ''
    # print(matrixForConsensusInteractions)
    # Sorting unitPairNames
    for loop_id in matrixForConsensusInteractions:
        for pairtype in matrixForConsensusInteractions[loop_id]:
            matrixForConsensusInteractions[loop_id][pairtype] = {key: matrixForConsensusInteractions[loop_id][pairtype][key] for key in sorted(matrixForConsensusInteractions[loop_id][pairtype])}

    # print(matrixForConsensusInteractions)
    return matrixForConsensusInteractions

def get_nts_by_columns(core_positions):
    transposed_nts = [list(x) for x in zip(*core_positions.values())]

    return transposed_nts

def motif_to_varna(core_positions = None, matrixForConsensusInteractions = None, name = None, output_dir = '.'):

    transposed_nts = [list(x) for x in zip(*core_positions.values())] # Get nts by columns
    # get sequenceDBN
    sequenceDBN = ''
    for each_column_nts in transposed_nts:
        # print(each_column_nts)
        each_column_nts = [item.split('|')[3] for item in each_column_nts]
        sequenceDBN = sequenceDBN + analyze_sequence_column(each_column_nts)
    
    #get structureDBN
    # auxBPs, structureDBN = get_auxBPs_and_structureDBN(matrixForConsensusInteractions, len(sequenceDBN))

    # This version of varna will change to 3.7 when it goes to server
    # Attection! In server spaceBetweenBases is no longer needed
    # print(sequenceDBN, structureDBN)
    # varna_command = 'java -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -algorithm "naview" -spaceBetweenBases 1.4 -sequenceDBN "%s" -structureDBN "%s" -baseNum "#334455" -periodNum 1 -auxBPs "%s" -o clusters/varnaImages/%s.png' % (sequenceDBN, structureDBN, auxBPs, name)
    # print(varna_command)

    # use os.command instead
    # os.system(varna_command)
    auxBPs, structureDBN = get_most_common_interaction(matrixForConsensusInteractions, len(sequenceDBN))
    # print('get_most_common_interaction from matrix: %s' % matrixForConsensusInteractions['most_common_interactions'])
    # matrixForConsensusInteractions['sequenceDBN'] = sequenceDBN
    # sequenceDBN, structureDBN, auxBPs, name
    varna_command_info = {
        'sequenceDBN': sequenceDBN,
        'structureDBN': structureDBN,
        'auxBPs': auxBPs,
        'output_dir': output_dir
    }
    # print('varna_command: %s' % varna_command)
    
    return varna_command_info

def is_key_in_dict(dictionary, key):
    return key in dictionary

def get_most_common_interaction(matrixForConsensusInteractions, loop_length):
    interaction_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    categories_to_consider = ['basepair', 'pairsStacks', 'BPh', 'BR', 'others']

    true_interaction_dict = defaultdict(lambda: defaultdict())
    # print(matrixForConsensusInteractions)

    cww_pairs =[]

    auxBPs = ''

    # get the true interaction for every unit pair in all categories
    for loop_id, categories in matrixForConsensusInteractions.items():
        for category in categories_to_consider:
            current_category_dict = categories[category]
            for unitPairName, interaction in current_category_dict.items():
                if interaction!='':
                    interaction_counts[category][unitPairName][interaction] += 1
                if not interaction.startswith('n') and interaction!='':
                    true_interaction_dict[category][unitPairName] = interaction

    most_common_interactions = {}
    numOfInstance = len(matrixForConsensusInteractions)

    for category, unitPairName_interactions in interaction_counts.items():
        most_common_interactions[category] = {}

        for unitPairName, interaction_and_counts in unitPairName_interactions.items():
            most_common_interaction = ''
            true_interaction_count = 0
            near_interaction_count = 0
            # print('now! the category is %s \nThe unitPairName is %s' % (category, unitPairName))
            # print(true_interaction_dict)
            # print(interaction_counts)
            if is_key_in_dict(true_interaction_dict[category], unitPairName):
                for interaction, count in interaction_and_counts.items():
                    
                        if interaction == true_interaction_dict[category][unitPairName]:
                            true_interaction_count = count
                        elif interaction == 'n' + true_interaction_dict[category][unitPairName]:
                            near_interaction_count = count
                # criterion
                criterion_1 = true_interaction_dict[category][unitPairName]
                criterion_2 = (6*true_interaction_count + 4*near_interaction_count) > 3*numOfInstance
                # old criterion from matlab code
                # (3*T + 2*NC + NN > L) and (T + NC > L/8) or (3*T + 2*NC > 30)


                if criterion_1 and criterion_2:
                    most_common_interaction = true_interaction_dict[category][unitPairName]
                    if category == 'basepair':
                        # get auxBPs 
                        x, y = [int(i) - 1 for i in unitPairName.split(' - ')]
                        if most_common_interaction == 'cWW':
                            cww_pairs.append((x, y))
                        edge_mapping = {'W': 'wc', 'H': 'h', 'S': 's'}
                        stericity_mapping = {'c': 'cis', 't': 'trans'}
                        # print('unitPairName is %s, and most_common_bptype is %s' % (unitPairName, most_common_interaction))
                        # auxBPs += '{} edge5={},'.format(unitPairName.replace("-",":"), edge_mapping[most_common_interaction[1]])
                        auxBPs += '(%s):edge5=%s,' % (unitPairName.replace(" - ",","), edge_mapping[most_common_interaction[1]])
                        # auxBPs += 'edge3={},'.format(edge_mapping[most_common_interaction[2]])
                        auxBPs += 'edge3=%s,' % (edge_mapping[most_common_interaction[2]])
                        # auxBPs += 'stericity={}; '.format(stericity_mapping[most_common_interaction[0]])
                        auxBPs += 'stericity=%s;' % (stericity_mapping[most_common_interaction[0]])
            
            most_common_interactions[category][unitPairName] = most_common_interaction
    matrixForConsensusInteractions['most_common_interactions'] = most_common_interactions
    # Generate structureDBN
    structureDBN = ['.' for _ in range(loop_length)]
    # print("loop_length: %s" % loop_length)
    for x, y in cww_pairs:
        structureDBN[x] = '('
        structureDBN[y] = ')'
    structureDBN = ''.join(structureDBN)
    # print(auxBPs)

    return auxBPs.rstrip(), structureDBN   
        

def find_bp_signature(queries, core_positions, centroid_id, matrixForConsensusInteractions):

    most_common_interactions = matrixForConsensusInteractions['most_common_interactions']['basepair']

    strand_info = queries[centroid_id]["loop_info"]["strand"]
    id_and_its_strand = []
    for unit_id in core_positions[centroid_id]:
        for strand in strand_info:
            if unit_id in strand:
                id_and_its_strand.append([unit_id, strand_info.index(strand)])



    dim_bp_matrix = len(core_positions[centroid_id])

    bp_matrix = [[0 for _ in range(dim_bp_matrix)] for _ in range(dim_bp_matrix)]

    paired_units = set()
    not_null_interactions = set()
    for unit_pair in most_common_interactions:
        if most_common_interactions[unit_pair] != '':
            x, y = [int(i) for i in unit_pair.split(' - ')]
            x, y = x-1, y-1
            bp_matrix[x][y] = most_common_interactions[unit_pair]
            paired_units.add(x)
            paired_units.add(y)
            not_null_interactions.add(unit_pair)

    unpaired_units = set()
    for i in range(dim_bp_matrix):
        if i not in paired_units:
            unpaired_units.add(i)

    bpSignature = []

    try:
        first_cww = '1 - %s' %(dim_bp_matrix)
        if most_common_interactions[first_cww] == 'cWW':
            bpSignature.append('cWW')
            not_null_interactions.remove(first_cww)
    except:
        print('Debugging info: An error happened on %s' % centroid_id)


    # scan the matrix
    for row in range(dim_bp_matrix):
        # scan by row
        for unit in range(dim_bp_matrix):
            if bp_matrix[row][unit] != 0:
                unit_pair = '%s - %s' % (row + 1, unit + 1)
                if unit_pair in not_null_interactions:
                    bpSignature.append(bp_matrix[row][unit])
                    not_null_interactions.remove(unit_pair)

            if row in unpaired_units:
                bpSignature.append("L")
                unpaired_units.remove(row)

        
        
        # scan by column
        column = range(dim_bp_matrix)[-row - 1]
        for another_row in reversed(range(dim_bp_matrix)):
            if bp_matrix[another_row][column] != 0:
                unit_pair = '%s - %s' % (another_row + 1, column + 1)
                if unit_pair in not_null_interactions:
                    bpSignature.append(bp_matrix[another_row][column])
                    not_null_interactions.remove(unit_pair)

            if column in unpaired_units:
                bpSignature.append("R")
                unpaired_units.remove(column)
                
    if bpSignature[0] != 'cWW':
        try:
            switch = bpSignature.index('cWW')
            bpSignature[0], bpSignature[switch] = 'cWW', bpSignature[0]
        except:
            pass
    # should be not startwith IL, then it will work for J3 as well 
    if not centroid_id.startswith('IL_'):
        for i in bpSignature:
            if i == "L" or i == "R":
                bpSignature[bpSignature.index(i)] = "F"
    bpSignature = "-".join(bpSignature)

    return bpSignature





def valid_list(lst):
    return all(elem in ['A', 'C', 'G', 'U'] for elem in lst)


def analyze_sequence_column(nts):
    #ignore them if not acgu  dont crash 
    if not valid_list(nts):
        print("This posistion has some other nucleotides besides A, C, G, U")
        # raise ValueError('Unexpected case: nts list is wrong')

    A = 'A' in nts
    C = 'C' in nts
    G = 'G' in nts
    U = 'U' in nts

    if not A and not G and C and U:
        consensus = 'Y'  # pyrimidines
    elif A and G and not C and not U:
        consensus = 'R'  # purines
    elif A and G and C and U:
        consensus = 'N'  # any nucleotide
    elif A and not G and not C and not U:
        consensus = 'A'
    elif not A and G and not C and not U:
        consensus = 'G'
    elif not A and not G and C and not U:
        consensus = 'C'
    elif not A and not G and not C and U:
        consensus = 'U'
    elif A and not G and not C and U:
        consensus = 'W'  # weak
    elif not A and G and C and not U:
        consensus = 'S'  # strong
    elif not A and G and not C and U:
        consensus = 'K'  # keto
    elif A and not G and C and not U:
        consensus = 'M'  # amino
    elif A and G and not C and U:
        consensus = 'D'  # not C
    elif A and G and C and not U:
        consensus = 'V'  # not U
    elif A and not G and C and U:
        consensus = 'H'  # not G
    elif not A and G and C and U:
        consensus = 'B'  # not A
    else:
        raise ValueError('Unexpected case: cannot find a consensus')

    return consensus





# java -Djava.awt.headless=true -cp VARNAv3-7.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN "URYR" -structureDBN "(..)" -baseNum "#ffffff" -auxBPs "(1,4):edge5=wc,edge3=wc,stericity=cis;(2,3):edge5=wc,edge3=wc,stericity=cis;" -o MotifAtlas/Releases/ilmarch_redo/Groups/html/2ds/Group_169.png


def main():
    # fr3d_interaction_to_pair_list = fr3d_interactions.get_fr3d_interaction_to_pair_list('4TNA')
    # lst_all_types = set()
    # all_stacks= ['s35','s55','s33']
    # all_bptypes= ['cWW', 'tWW', 'cWH', 'tWH', 'cWS', 'tWS', 'cHH', 'tHH', 'cHS', 'tHS', 'cSS', 'tSS']
    # for key in fr3d_interaction_to_pair_list:
    #     if key not in (all_stacks + all_bptypes):
    #         lst_all_types.add(key)
    
    # print(lst_all_types)
    fr3d_pair_to_interaction_list, fr3d_pair_to_crossing_number = fr3d_interactions.get_fr3d_pair_to_interaction_list('4TNA')

    print("interaction_list:")
    print(fr3d_pair_to_interaction_list)
    print("crossing_number:")
    print(fr3d_pair_to_crossing_number)


if __name__ == "__main__":
    res = main()
