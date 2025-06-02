

from collections import defaultdict
import sys

# Modified nucleotide mappings
from fr3d.modified.mapping import modified_base_to_parent

if sys.version_info[0] < 3:
    import urllib
else:
    import urllib.request


def get_pdb_id(loop_id):
    pdb_id = loop_id.split('_')[1]
    return pdb_id


def get_nts_by_columns(core_positions):
    transposed_nts = [list(x) for x in zip(*core_positions.values())]

    return transposed_nts


def get_parent_as_RNA(sequence):
    """
    Look up parent sequence for RNA, DNA, and modified nucleotides.
    Return A, C, G, U to simplify getting a consensus.
    """

    if sequence in ['A','C','G','U']:
        return sequence
    elif sequence in ['DA','DC','DG']:
        return sequence[1]
    elif sequence == 'DT':
        return 'U'
    elif sequence in modified_base_to_parent:
        parent = modified_base_to_parent[sequence]
        if parent in ['A','C','G','U']:
            return parent
        elif parent == 'DT':
            return 'U'
        elif parent in ['DA','DC','DG']:
            return parent[1]

    return ''


def analyze_sequence_column(nts):
    """
    nts is a list of parent nucleotides like A, C, G, U
    """

    if len(nts) == 0:
        return 'N'

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
        # avoid crashing in every case
        consensus = '?'

    return consensus


def get_consensus_interactions(matrixForConsensusInteractions, loop_length):
    """
    criterion1 and criterion2 are for consensus interactions
    Also returns VARNA commands
    """

    interaction_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    categories_to_consider = ['basepair', 'pairsStacks', 'BPh', 'BR', 'others']

    category_pair_to_true_interaction = defaultdict(lambda: defaultdict())

    cww_pairs = []

    auxBPs = ''

    # get the true interaction for every position pair in all categories
    for loop_id, categories in matrixForConsensusInteractions.items():
        for category in categories_to_consider:
            current_category_dict = categories[category]
            for corePositionPair, interaction in current_category_dict.items():
                if interaction != '':
                    interaction_counts[category][corePositionPair][interaction] += 1
                if not interaction.startswith('n') and interaction != '':
                    category_pair_to_true_interaction[category][corePositionPair] = interaction

    consensus_interactions = {}
    num_instances = len(matrixForConsensusInteractions)

    for category, corePositionPair_interactions in interaction_counts.items():
        consensus_interactions[category] = {}

        for corePositionPair, interaction_and_counts in corePositionPair_interactions.items():
            consensus = ''
            true_interaction_count = 0
            near_interaction_count = 0

            if corePositionPair in category_pair_to_true_interaction[category]:
                for interaction, count in interaction_and_counts.items():

                    if interaction == category_pair_to_true_interaction[category][corePositionPair]:
                        true_interaction_count = count
                    elif interaction == 'n' + category_pair_to_true_interaction[category][corePositionPair]:
                        near_interaction_count = count

                # criteria to decide if the interaction occurs enough times to be conserved
                criterion_1 = category_pair_to_true_interaction[category][corePositionPair]
                criterion_2 = (6*true_interaction_count + 4*near_interaction_count) > 3*num_instances
                # old criterion from matlab code
                # (3*T + 2*NC + NN > L) and (T + NC > L/8) or (3*T + 2*NC > 30)

                if criterion_1 and criterion_2:
                    consensus = category_pair_to_true_interaction[category][corePositionPair].replace("n","")
                    if category == 'basepair':
                        # get auxBPs for VARNA
                        x, y = [int(i) - 1 for i in corePositionPair.split(' - ')]
                        if consensus == 'cWW':
                            cww_pairs.append((x, y))
                        edge_mapping = {'W': 'wc', 'H': 'h', 'S': 's'}
                        stericity_mapping = {'c': 'cis', 't': 'trans'}

                        edge1 = consensus[1].upper()
                        if not edge1 in edge_mapping:
                            continue
                        edge2 = consensus[2].upper()
                        if not edge2 in edge_mapping:
                            continue
                        stericity = consensus[0]
                        if not stericity in stericity_mapping:
                            continue
                        # print('corePositionPair is %s, and most_common_bptype is %s' % (corePositionPair, consensus))
                        # auxBPs += '{} edge5={},'.format(corePositionPair.replace("-",":"), edge_mapping[consensus[1]])
                        auxBPs += '(%s):edge5=%s,' % (corePositionPair.replace(" - ",","), edge_mapping[edge1])
                        # auxBPs += 'edge3={},'.format(edge_mapping[consensus[2]])
                        auxBPs += 'edge3=%s,' % (edge_mapping[edge2])
                        # auxBPs += 'stericity={}; '.format(stericity_mapping[consensus[0]])
                        auxBPs += 'stericity=%s;' % (stericity_mapping[stericity])

            consensus_interactions[category][corePositionPair] = consensus

    # Generate structureDBN ... DBN = dot-bracket notation
    structureDBN = ['.' for _ in range(loop_length)]

    # force the loop to start and end with a WC basepair, even if near
    # if not (1,loop_length) in cww_pairs:
    #     cww_pairs.append((1,loop_length))

    for x, y in cww_pairs:
        structureDBN[x] = '('
        structureDBN[y] = ')'
    # force the loop to start and end with a WC basepair, even if near
    # structureDBN[0] = '('
    # structureDBN[-1] = ')'
    # how to insist that the end of an IL or J strand pairs with the next one?  Need )( in the structureDBN

    structureDBN = ''.join(structureDBN)

    return auxBPs.rstrip(), structureDBN, consensus_interactions


def motif_to_varna(core_positions = None, matrixForConsensusInteractions = None, name = None, output_dir = '.'):

    transposed_nts = [list(x) for x in zip(*core_positions.values())] # Get nts by columns

    # get sequence summary letter by column
    sequenceDBN = ''
    for column_nts in transposed_nts:
        each_column_nts = []
        for item in column_nts:
            parent = get_parent_as_RNA(item.split('|')[3])
            if parent:
                each_column_nts.append(parent)
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
    auxBPs, structureDBN, consensus_interactions = get_consensus_interactions(matrixForConsensusInteractions, len(sequenceDBN))
    # matrixForConsensusInteractions['sequenceDBN'] = sequenceDBN
    # sequenceDBN, structureDBN, auxBPs, name
    varna_command_info = {
        'sequenceDBN': sequenceDBN,
        'structureDBN': structureDBN,
        'auxBPs': auxBPs,
        'output_dir': output_dir
    }
    # print('varna_command: %s' % varna_command)

    return varna_command_info, consensus_interactions


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

    # from pymotifs.motif_atlas.fr3d_interactions import get_fr3d_pair_to_interaction_list
    # fr3d_pair_to_interaction_list, fr3d_pair_to_crossing_number = get_fr3d_pair_to_interaction_list('4TNA')

    # print("interaction_list:")
    # print(fr3d_pair_to_interaction_list)
    # print("crossing_number:")
    # print(fr3d_pair_to_crossing_number)

    pass

if __name__ == "__main__":
    res = main()
