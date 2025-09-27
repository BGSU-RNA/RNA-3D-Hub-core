
from collections import defaultdict
import csv
import json
import networkx as nx
import numpy as np
import os
import os.path
import pickle
import sys
from sys import path
import scipy.io as sio
from time import time

# import the version of urlretrieve appropriate to the Python version
if sys.version_info[0] < 3:
    from urllib import urlretrieve as urlretrieve
else:
    from urllib.request import urlretrieve as urlretrieve

from pymotifs.motif_atlas.fr3d_interactions import get_fr3d_pair_to_interaction_list

DATAPATHATLAS      = '/usr/local/pipeline/hub-core/MotifAtlas/tier2/'
DATAPATHLOOPS      = '/usr/local/pipeline/hub-core/MotifAtlas/tier2/loops'
DATAPATHMOTIFDATA = os.path.join(DATAPATHATLAS,'motif_data')


# def download_and_save_motif_atlas_release_json_file(loop_type,release,path=DATAPATHATLAS):
#     """
#     Download a release of the motif atlas
#     Download json file locally in DATAPATHLOOPS
#     Save as pickle file
#     """

#     if not os.path.exists(path):
#         os.mkdir(path)

#     local_file_name = "%s_motif_atlas_release_%s.json" % (loop_type,release)
#     local_file_path = os.path.join(path,local_file_name)
#     # url = "https://rna.bgsu.edu/rna3dhub/motifs/release/%s/%s/json" % (loop_type,release)

#     # print("Downloading motif atlas release from %s to %s" % (url,path))
#     # urlretrieve(url,local_file_path)
#     print("Downloading motif atlas release %s to %s" % (local_file_name,path))
#     url = "https://rna.bgsu.edu/rna3dhub/motifs/release/%s/%s/json" % (loop_type,release)
#     print(url)
#     print(local_file_path)
#     try:
#         urlretrieve(url,local_file_path)
#     except IOError:
#         print('Cannot download motif atlas json file %s from production' % url)

#     loops_and_strands = json.load(open(local_file_path, "r"))

#     pickle_file_name = "%s_motif_atlas_release_%s.pickle" % (loop_type,release)
#     pickle_file_path = os.path.join(path,pickle_file_name)
#     pickle.dump(obj=loops_and_strands,file=open(pickle_file_path,"wb"),protocol=2)

#     return

# def load_motif_atlas_release(loop_type,release='3.57',path=DATAPATHATLAS):
#     """
#     Download a release of the motif atlas
#     """
#     file_name = loop_type + "_motif_atlas_release_" + release + ".pickle"
#     file_path = os.path.join(path,file_name)

#     if not os.path.exists(file_path):
#         download_and_save_motif_atlas_release_json_file(loop_type,release,path)

#     loops_and_strands = []
#     loops_and_strands = pickle.load(open(file_path,"rb"))

#     return loops_and_strands


def get_highest_and_median_disc(search_results,loops):
    if len(loops)<=2:
        return 0,0
    highest_disc = 0
    median_disc = 0
    disc_list = []
    for i in range(len(loops)):
        for j in range(len(loops)):
            if i == j:
                continue
            loops_pair = (loops[i],loops[j])
            if loops_pair in search_results.keys():
                if search_results[loops_pair][0]["discrepancy"] < 99:
                    disc_list.append(search_results[loops_pair][0]["discrepancy"])

    highest_disc = max(disc_list)

    disc_list.sort()
    mid = len(disc_list) // 2
    median_disc = (disc_list[mid] + disc_list[~mid]) / 2
    return(highest_disc,median_disc)

def get_dq_between_motif_groups(cliques,dq_matrix,groups):
    """
    Testing function, used in cluster_loops
    """
    cliques.sort(key=len,reverse=True)
    dq_between_motifs ={}
    for clique_i in range(len(cliques)):
        for clique_j in range(clique_i+1,len(cliques)):
            dq_codes = set()
            for i in cliques[clique_i]:
                for j in cliques[clique_j]:
                    dq_codes.add(dq_matrix[i][j])
            dq_between_motifs[(clique_i,clique_j)]=dq_codes
    dq_between_motifs
    for motifs_pair in dq_between_motifs.keys():
        print("The disqualification codes between %s is %s" % (motifs_pair,dq_matrix[motifs_pair]))
    return(dq_between_motifs)

def print_discrepancy_matrix(matrix,all_loops_ids):
    """'
    Helper function
    """
    df = pd.DataFrame(matrix,index = all_loops_ids)
    with pd.option_context('display.max_rows', None,
                       'display.max_columns', None,
                       'display.precision', 3,
                       ):
        # print(df)
        pass
    return()

def show_discrepancy(search_results,targets,loops_of_interest):
    dq_codes = defaultdict(list)
    discs = []
    highest_disc = 0
    highest_disc_loop = ""
    for target in targets:
        for loop in loops_of_interest:
            if loop == target:
                continue
            if (target,loop) in search_results.keys():
                dq = search_results[(target,loop)][0]['dq']
                disc = search_results[(target,loop)][0]['discrepancy']
                if dq == [-1] or dq ==[0]:
                    dq = search_results[(loop,target)][0]['dq']
                    disc = search_results[(loop,target)][0]['discrepancy']
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
        print("Highest disc:{}, Average disc:{}".format(highest_disc,average_disc))

def validate_clusters(loop_type,cluster_method,search_results,loops_of_interest,python_clustered_motif_groups,path = "./data/"):
    """
    Calculate statistics and write to output file
    """

    total_identical_groups = 0
    file_name = "%s_%s_validation.txt" % (loop_type,cluster_method)
    file_path = os.path.join(path,file_name)
    file_path = os.path.join(file_name)         # save locally, easier to find
    file1 = open(file_path,"w")
    motif_url = "https://rna.bgsu.edu/rna3dhub/motif/view/"
    motif_atlas_release = load_motif_atlas_release(loop_type)
    matlab_motif_to_loops,loops_to_matlab_motif = motif_ids(motif_atlas_release,loops_of_interest)
    if len(matlab_motif_to_loops) != len(python_clustered_motif_groups):
        file1.write("The amount of groups are incorrect")
    for matlab_motif,matlab_group in sorted(matlab_motif_to_loops.items(),key=lambda k: len(k[1])):
        new_motif = True #Bool to print
        matlab_group.sort()
        for python_group in python_clustered_motif_groups:
            python_group.sort()
            common_loops = set(python_group) & set(matlab_group)
            if python_group == matlab_group:
                total_identical_groups +=1
                file1.write("Motif {} with length: {} agree with Python".format(matlab_motif,len(matlab_group)))
                if len(python_group) >1:
                    file1.writelines(python_group)
            if len(common_loops)> 0  and python_group != matlab_group:
                if new_motif:
                    highest_disc,median_disc = get_highest_and_median_disc(search_results,matlab_group)
                    file1.write("=========================================================================\n")
                    file1.write("Motif group: {}. Len : {}. Highest disc: {}. Median disc: {}. {}\n".format(matlab_motif,len(matlab_group),highest_disc,median_disc,motif_url+matlab_motif))
                    file1.write("\n")
                    file1.write("Matlab: {}\n".format(matlab_group))
                    new_motif = False
                    file1.write("\n")
                    for i in range(len(matlab_group)):
                        for j in range(len(matlab_group)):
                            if i ==j or (matlab_group[i],matlab_group[j]) not in search_results.keys():
                                continue
                            dq = search_results[(matlab_group[i],matlab_group[j])][0]['dq']
                            if len(dq) > 0 and -1 not in dq and 0 not in dq:
                                file1.write("{} - {}: {}".format(matlab_group[i],matlab_group[j],dq))
                file1.write("\n")
                highest_disc,median_disc = get_highest_and_median_disc(search_results,python_group)
                file1.write("Python. Highest disc: {}. Median disc: {}. {}\n".format(highest_disc,median_disc,python_group))
                # for i in range(len(python_group)):
                #         for j in range(len(python_group)):
                #             if i ==j or (python_group[i],python_group[j]) not in search_results.keys():
                #                 continue
                #             dq = search_results[(python_group[i],python_group[j])][0]['dq']
                #             if len(dq) > 0 and -1 not in dq and 0 not in dq:
                #                 print("{} - {}: {}".format(python_group[i],python_group[j],dq))
                file1.write(" \n")
                file1.write("\tIntersection: {}\n".format(set(python_group) & set(matlab_group)))
                file1.write("\tLoops in Python but not Matlab:")
                differences = defaultdict(list)
                for loop in (set(python_group)-set(matlab_group)):
                    if(loop in loops_to_matlab_motif.keys()):
                        differences[loops_to_matlab_motif[loop]].append(loop)
                for motif,loops in sorted(differences.items(),key = lambda k: k[0]):
                    file1.write("\tFrom motif {} - {}".format(motif,motif_url+motif))
                    file1.write("\t{}".format(sorted(loops)))

                file1.write("\n")

                file1.write("\tLoops in Matlab but not Python:")
                differences = defaultdict(list)
                for loop in (set(matlab_group)-set(python_group)):
                    differences[loops_to_matlab_motif[loop]].append(loop)
                for motif,loops in sorted(differences.items(),key = lambda k: k[0]):
                    file1.write("\tFrom motif {} - {}".format(motif,motif_url+motif))
                    file1.write("\t{}".format(sorted(loops)))
                file1.write("\n")
    file1.write("Total identical groups: {}".format(total_identical_groups))
    return()

def motif_ids(motif_atlas_release, loops_of_interest = False):
    """
    Map loop id to motif id, and motif id to list of loop ids
    Return a dictionary: {loop_id:motif_id}
    For testing.
    """
    motif_id_to_loop_ids = defaultdict(list)
    loop_id_to_motif_id = defaultdict(list)

    if(loops_of_interest):
        for motif_group in motif_atlas_release:
            motif_id = motif_group["motif_id"]
            loop_alignments = motif_group["alignment"]
            if loop_alignments:
                for loop_id,alignment in loop_alignments.items():
                    if loop_id in loops_of_interest:
                        motif_id_to_loop_ids[motif_id].append(loop_id)
                        loop_id_to_motif_id[loop_id] = motif_id
    else:
        for motif_group in motif_atlas_release:
            motif_id = motif_group["motif_id"]
            loop_alignments = motif_group["alignment"]
            if loop_alignments:
                for loop_id,alignment in loop_alignments.items():
                    motif_id_to_loop_ids[motif_id].append(loop_id)
                    loop_id_to_motif_id[loop_id] = motif_id

    return motif_id_to_loop_ids,loop_id_to_motif_id

def heat_map(cliques_lst,groups,disc_m):
    """
    Testing function, used for cluster_loops.
    """
    best_cliques = []
    while True:
        cliques_lst.sort(key=len,reverse=True)
        if len(cliques_lst[0]) == 1 or len(cliques_lst) == 0 or len(cliques_lst[0]) == 0:
            break
        best_clique = find_best_clique(disc_m,cliques_lst,1)
        best_cliques.append(best_clique)
        # remove row/column from both matrices
        cliques_lst = [[x for x in sub if x not in best_clique] for sub in cliques_lst]
    for clique in range(len(best_cliques)):
        group_disc_matrix = np.zeros(shape=(len(best_cliques[clique]),len(best_cliques[clique])))
        for i in range(len(best_cliques[clique])):
            for j in range(i+1,len((best_cliques[clique]))):
                group_disc_matrix[i][j] = disc_m[best_cliques[clique][i]][best_cliques[clique][j]]
                group_disc_matrix[j][i] = disc_m[best_cliques[clique][j]][best_cliques[clique][i]]
        df = pd.DataFrame(group_disc_matrix,index=groups[clique],columns=groups[clique])
        sns.set(font_scale=0.5)
        sns.heatmap(df,annot=True)
        plt.show()
    return()

def load_all_clusters(loop_type,path = "./data/"):
    all_clusters = []
    file_name = loop_type + "_all_clusters.pickle"
    file_path = os.path.join(path,file_name)

    if os.path.exists(file_path):
        if sys.version_info[0] < 3:
            [search_results,all_loops_ids] = pickle.load(open(file_path, "rb"))
        else:
            all_clusters = pickle.load(open(file_path, "rb"), encoding = 'latin1')

    return(all_clusters)

def load_all_queries(loop_type,path = "./data/"):
    queries = []
    file_name = loop_type + "_all_queries.pickle"
    file_path = os.path.join(path,file_name)

    if os.path.exists(file_path):
        if sys.version_info[0] < 3:
            [search_results,all_loops_ids] = pickle.load(open(file_path, "rb"))
        else:
            queries = pickle.load(open(file_path, "rb"), encoding = 'latin1')

    else:
        print("Not able to load %s" % file_path)

    return(queries)

def load_dist_matrices(loop_type,path = "./data/"):
    file_name = loop_type + "_dist_matrices.pickle"
    file_path = os.path.join(path,file_name)

    if os.path.exists(file_path):
        if sys.version_info[0] < 3:
            [search_results,all_loops_ids] = pickle.load(open(file_path, "rb"))
        else:
            dist_matrices = pickle.load(open(file_path, "rb"), encoding = 'latin1')

    return(dist_matrices)

def save_all_clusters(loop_type,clusters,path = "./data/"):
    file_name = loop_type + "_all_clusters.pickle"
    file_path = os.path.join(path,file_name)
    if not os.path.exists(path):
        os.mkdir(path)
    pickle.dump(obj = clusters, file = open(file_path, "wb"), protocol = 2)
    return()

def save_all_queries(loop_type,queries,path = "./data/"):
    file_name = loop_type + "_all_queries.pickle"
    file_path = os.path.join(path,file_name)
    if not os.path.exists(path):
        os.mkdir(path)
    pickle.dump(obj = queries, file = open(file_path, "wb"), protocol = 2)
    return()

def save_all_search_results(loop_type,search_results,all_loops_ids,path = "./data/"):
    file_name = loop_type + "_all_search_results.pickle"
    file_path = os.path.join(path,file_name)
    if not os.path.exists(path):
        os.mkdir(path)
    pickle.dump(obj = [search_results,all_loops_ids], file = open(file_path, "wb"), protocol = 2)
    return()

def save_dist_matrices(loop_type,dist_matrices,path = "./data/"):
    file_name = loop_type + "_dist_matrices.pickle"
    file_path = os.path.join(path,file_name)
    if not os.path.exists(path):
        os.mkdir(path)
    pickle.dump(obj = dist_matrices, file = open(file_path, "wb"), protocol = 2)
    return()

def load_all_search_results(loop_type,path = "./data/"):
    search_results = {}
    all_loops_ids = []
    file_name = loop_type + "_all_search_results.pickle"
    file_path = os.path.join(path,file_name)

    if os.path.exists(file_path):
        if sys.version_info[0] < 3:
            [search_results,all_loops_ids] = pickle.load(open(file_path, "rb"))
        else:
            search_results,all_loops_ids = pickle.load(open(file_path, "rb"), encoding = 'latin1')

        for (query,search_space),value in search_results.items(): #Load No-match result as [{'dq': 0, 'discrepancy': 99}] instead of an integer
            if isinstance(value,int):
                search_results[(query,search_space)] = [{'dq': [value], 'discrepancy': 99}]
    return(search_results,all_loops_ids)

def findMainMatlabMotif(loop_type,candidates,loops_of_interest):
    """

    """
    sorted_candidates = sorted(candidates)

    motif_atlas_release = load_motif_atlas_release(loop_type)
    matlab_motif_to_loops,loops_to_matlab_motif = motif_ids(motif_atlas_release,loops_of_interest)
    differences= []
    motif_to_differences = {}
    max_common_loop_len = 0

    for matlab_motif,matlab_group in sorted(matlab_motif_to_loops.items(),key=lambda k: len(k[1]),reverse=True):
        matlab_group.sort()
        common_loops= set(sorted_candidates)&set(matlab_group)
        if sorted_candidates == matlab_group:
            return (motif_to_differences)
        if len(common_loops)> 0  and sorted_candidates != matlab_group:
            differences = list(set(matlab_group)-set(sorted_candidates))
            motif_to_differences[matlab_motif] = differences
    return(motif_to_differences)



"""
Making HTML output
"""


def writeHTMLOutput(num,loops_of_interest,cluster_method,search_results,motif_groups,loop_id_to_core_units,candidates,queries,loop_ids_to_discrepancy,loop_id_to_matching_groups,loop_id_to_data,loop_id_to_group_number,group_max_distance,release, loop_id_type_pair_to_interaction, varna_command_info, bpSignature, consensus_interactions):
    """
    Write the list of candidates in an HTML format that also shows
    the coordinate window and a heat map of all-against-all distances.
    Include information about matches to Matlab motif groups.
    num is the group number

    candidates is the list of loop ids in the current group

    """

    JS1 = '  <script src="./js/JSmol.min.nojq.js"></script>'
    JS2 = '  <script src="./js/jquery.jmolTools.js"></script>'
    JS3 = '  <script src="./js/imagehandlinglocal.js"></script>'
    JS4 = '<script src="./js/jmolplugin.js" type="text/javascript"></script>'
    JS5 = '<script type="text/javascript" src="./js/heatmap.js"></script>'
    TEMPLATEPATH = './search/'
    OUTPUTPATH = "clusters/"

    print("Writing HTML file for group %3d with %3d candidates" % (num,len(candidates)))

    Q = defaultdict(list)
    Q["name"] = "%03d group" % num
    Q["numFilesSearched"] = 5
    Q["searchFiles"]= 'searchFiles'
    Q["elapsedCPUTime"] = 5

    loop_type = loops_of_interest[0].split("_")[0]

    motif_atlas_release = load_motif_atlas_release(loop_type,release)
    matlab_motif_to_loops,loops_to_matlab_motif = motif_ids(motif_atlas_release,loops_of_interest)
    motif_to_differences = findMainMatlabMotif(loop_type, candidates, loops_of_interest)

    pagetitle = "%s" % Q['name']
    motif_url = "https://rna.bgsu.edu/rna3dhub/motif/view/"
    #htmlfilename = Q['name'].replace(" ","_")

    # Create links to previous and next motif groups
    buttons = '<div class="btn-group">\n'
    if num == 0:
        buttons += '<a href = "%s_%s_%03d.html">\n' % (loop_type,cluster_method,num)
    else:
        buttons += '<a href = "%s_%s_%03d.html">\n' % (loop_type,cluster_method,num-1)
    buttons += 'Previous</a> | \n'
    if num ==len(motif_groups)-1:
        buttons += '<a href = "%s_%s_%03d.html">\n' % (loop_type,cluster_method,num)
    else:
        buttons += '<a href = "%s_%s_%03d.html">\n' % (loop_type,cluster_method,num+1)
    buttons += 'Next</a>\n'
    buttons += "</div>\n"
    buttons += "Maximum discrepancy within this group is %8.4f" % group_max_distance[num]

    # build the table that lists the candidates
    candidatelist = '<table style="white-space:nowrap;">\n'

    # write header line
    candidatelist += "<tr><th>S.</th><th>Show</th>"

    candidatelist += "<th>Loop ID</th>"
    candidatelist += "<th>Matlab group</th>"
    candidatelist += "<th>Annotation</th>"
    candidatelist += "<th># of core nts</th>"
    candidatelist += "<th># of non-core nts</th>"
    candidatelist += "<th># of bulged nts</th>"
    candidatelist += "<th>Could be in</th>"

    # for i in range(0,len(loop_id_to_core_units[candidates[0]])):
    #     candidatelist += "<th></th>"
    #     candidatelist += "<th>" + str(i+1)+"</th>"
    # use the largest number of loop to write header
    largestNumOfLoop = max(len(loop_id_to_core_units[loop_id]) for loop_id in loop_id_to_core_units)
    for i in range(0,largestNumOfLoop):
        candidatelist += "<th></th>"
        candidatelist += "<th>" + str(i+1)+"</th>"


    # add consensus interactions
    pairsToPrint = defaultdict(list)

    pairTypes = ['basepair','pairsStacks','BPh','BR','others']

    for pairtype in pairTypes:
        for corePositionPair, consensus_interaction in loop_id_type_pair_to_interaction[candidates[0]][pairtype].items():
            candidatelist += "<th>%s</th>" % (corePositionPair)

    candidatelist += "</tr>\n"

    # for corePositionPair, consensus_interaction in loop_id_type_pair_to_interaction[candidates[0]].items():
    #     candidatelist += "<td>%s</td>" % (corePositionPair)

    # write one row for each candidate
    for i in range(0,len(candidates)):
        candidate = candidates[i]
        highest_disc,mean_disc = loop_id_to_matching_groups[candidate][num]

        unitIDs= loop_id_to_core_units[candidate]
        #core_nts = sum(len(nts) for nts in queries[candidate]["loop_info"]["strand"])
        lengthOfLoop = sum(len(nts) for nts in queries[candidate]["loop_info"]["strand"])
        bulged_nts = len(queries[candidate]["loop_info"]["bulged"])
        core_nts = len(unitIDs)
        non_core_nts = lengthOfLoop - core_nts
        candidatelist += '<tr><td>'+str(i+1)+'.</td><td><label><input type="checkbox" id="'+str(i)+'" class="jmolInline" data-coord="'
        for j in range(0,core_nts):
            candidatelist += unitIDs[j] # append all unit ids for the current candidate
            if j < core_nts-1:
                candidatelist += ','
        candidatelist += '">&nbsp</label></td>'
        candidatelist += '<td><a href="https://rna.bgsu.edu/rna3dhub/loops/view/%s" target="_blank" rel="noopener noreferrer"> %s </a></td>' % (candidate,candidate)
        if not loops_to_matlab_motif[candidate]:
            loops_to_matlab_motif[candidate]=""
        url = motif_url+loops_to_matlab_motif[candidate]
        candidatelist += '<td><a href="%s" target="_blank" rel="noopener noreferrer">%s</a> (%.4f,%.4f)</td>' % (url,loops_to_matlab_motif[candidate],highest_disc,mean_disc)
        candidatelist += "<td>%s</td>" % (loop_id_to_data[candidate])
        candidatelist += "<td>%s</td>" % (core_nts)
        candidatelist += "<td>%s</td>" % (non_core_nts)
        candidatelist += "<td>%s</td>" % (bulged_nts)
        candidatelist += "<td>"

        # what other group could this loop be in?
        # why does this only ever list one???
        for group_num in sorted(loop_id_to_matching_groups[candidate].keys()):
            highest_disc, mean_disc = loop_id_to_matching_groups[candidate][group_num]
            candidatelist += '<a href = "%s_%s_%03d.html" target="_blank"> %03d' % (loop_type,cluster_method,group_num,group_num)
            if num == group_num:
                candidatelist += '</a> [%0.4f,%0.4f]' % (highest_disc,mean_disc)   # numbers for the current group
            else:
                candidatelist += '</a> (%0.4f,%0.4f)' % (highest_disc,mean_disc)
            candidatelist += '</br>'
        candidatelist += '</td>'

        for j in range(0, largestNumOfLoop):
            # print("largest num of loop: %s \n core_nts: %s \n j: %s" % (largestNumOfLoop, core_nts, j))
            # for j in range(0,core_nts):
            if j < core_nts:
                # print("j: %s" % j)
                unitID = unitIDs[j]
                unit_order = unitID.split("|")[4]
                unit_nt = unitID.split("|")[3]
                candidatelist += "<td>%s</td>" % (unit_nt)
                candidatelist += "<td>%s</td>" % (unit_order)
            else:
                candidatelist += "<td></td>"
                candidatelist += "<td></td>"

        pairTypes = ['basepair','pairsStacks','BPh','BR','others']
        for pairtype in pairTypes:
            for corePositionPair, consensus_interaction in loop_id_type_pair_to_interaction[candidate][pairtype].items():
                if consensus_interaction == '':
                    candidatelist += "<td></td>"
                else:
                    candidatelist += "<td>%s</td>" % (consensus_interaction)

        candidatelist += '</tr>\n'
    candidatelist += '<tr>\n'
    for i in range(0,9):
        if i != 2:
            candidatelist += '<td></td>'
        else:
            candidatelist += '<td>Consensus</td>'
    for i in range(len(varna_command_info['sequenceDBN'])):
        candidatelist += "<td>%s</td>" % (varna_command_info['sequenceDBN'][i])
        candidatelist += "<td></td>"

    for category, corePositionPair_interactions in consensus_interactions.items():
        for corePositionPair, most_common_interaction in corePositionPair_interactions.items():
            # print(f"For {corePositionPair}, the most common interaction is {most_common_interaction}")
            candidatelist += "<td>%s</td>" % (most_common_interaction)

    candidatelist += '</tr>\n'
    candidatelist += '</table>\n'


    # write discrepancy data in new 2022 list format
    # first element is a reference to the div in which the heatmap should appear
    discrepancydata = 'var data = ["#chart",[\n'            # start a list, start a matrix

    # second element is a matrix with the numerical values of the discrepancy
    # writing both upper and lower triangles of the matrix
    s = len(candidates)
    for c in range(0,s):
        discrepancydata += '['     # start a row of the discrepancy matrix
        loop1 = candidates[c]       # label rows with loop id
        for d in range(0,s):
            loop2 = candidates[d]
            discrepancydata += "%.4f" % loop_ids_to_discrepancy[loop1][loop2]  # one entry
            if d < s-1:
                discrepancydata += ','  # commas between entries in a row
            else:
                discrepancydata += '],\n'  # end a row, newline

    discrepancydata += '],\n'           # end the matrix, continue the list

    # third element is a list of labels of instances
    discrepancydata += '['              # start list of instances
    for c in range(0,s):
        loop_id = candidates[c]
        discrepancydata += '"' + loop_id + '"'    # write one instance name in quotes
        if c < s-1:
            discrepancydata += ","  # commas between instances
        else:
            discrepancydata += "]\n]" # end list of instances, end list of data

    # read template.html into one string
    with open(TEMPLATEPATH + 'template.html', 'r') as myfile:
        template = myfile.read()

    # replace ###PAGETITLE### with pagetitle
    template = template.replace("###PAGETITLE###",pagetitle)

    sequence_column = 0   # column to set in fixed width font
    template = template.replace("###sequencecolumn###",str(sequence_column))

    queryNote = "Motif group %03d using %s" % (num,cluster_method)

    if "moreCandidatesThanHeatMap" in Q:
        queryNote += " " + Q["moreCandidatesThanHeatMap"] + "\n"
    else:
        queryNote += "\n"

    seeModifyQuery = ""
    template = template.replace("###QUERYNAME###",queryNote)
    for matlab_motif, differences in sorted(motif_to_differences.items(),key=lambda k: len(k[1]),reverse=True):
        if len(differences) > 0:
            seeModifyQuery += "<p> %d loops in Matlab motif <a href='%s'>%s</a> but not in this Python group: " % (len(differences),motif_url + matlab_motif,matlab_motif)
            for loop in differences:
                seeModifyQuery += '<a href="https://rna.bgsu.edu/rna3dhub/loops/view/%s" target="_blank">%s</a> ' % (loop,loop)
                if loop in loop_id_to_group_number:
                    gn = loop_id_to_group_number[loop]
                    seeModifyQuery += '(<a href="%s_%s_%03d.html" target="_blank">%03d</a>), ' % (loop_type,cluster_method,gn,gn)
                else:
                    seeModifyQuery += '(---), '

            seeModifyQuery += "</p>\n"

    template = template.replace("###SEEMODIFYQUERY###",seeModifyQuery)
    template = template.replace("###seeCSVOutput###",buttons)
    template = template.replace("###CANDIDATELIST###",candidatelist)
    template = template.replace("###JS1###",JS1)
    template = template.replace("###JS2###",JS2)
    template = template.replace("###JS3###",JS3)
    template = template.replace("###JS4###",JS4)

    # need to change
    varna_command = 'java -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -algorithm "naview" -spaceBetweenBases 1.4 -sequenceDBN "%s" -structureDBN "%s" -baseNum "#334455" -periodNum 1 -auxBPs "%s" -o clusters/varnaImages/%s_%s_%03d.png' % (varna_command_info['sequenceDBN'], varna_command_info['structureDBN'], varna_command_info['auxBPs'], loop_type,cluster_method,num)
    os.system(varna_command)

    template += '<img src="varnaImages/%s_%s_%03d.png" alt="varna image" style="height:600px;">' % (loop_type,cluster_method,num)
    # bpSignature = find_bp_signature(queries, loop_id_to_core_units, centroid_loop, loop_id_type_pair_to_interaction)
    template += "<p> basepair signature: %s </p>\n" % (bpSignature)



    refresh = ""
    if "reloadOutputPage" in Q and Q["reloadOutputPage"]:
        refresh = '<meta http-equiv="refresh" content="%d">' % REFRESHTIME
    template = template.replace("###REFRESH###",refresh)

    if len(candidates) > 1:
        template = template.replace("###JS5###",JS5)    # include heatmap.js code
        discrepancydata = '<script type="text/javascript">\n' + discrepancydata + '\n</script>'
        template = template.replace("###DISCREPANCYDATA###",discrepancydata)
    else:
        #template = template.replace("###DISCREPANCYDATA###","")
        template = template.replace("###JS5###","")    # do not display a heat map

    outputfilename = os.path.join(OUTPUTPATH,"%s_%s_%03d.html" % (loop_type,cluster_method,num))

    messages = ""

    messages += ",".join(candidates) + '<br>\n'

    messages += "\n<br>"
    if len(Q["userMessage"]) > 0:
        messages += "User messages:<br>\n"
        for line in Q["userMessage"]:
            messages += line + "<br>\n"
    else:
        messages += "No error or warning messages.<br>\n"

    template = template.replace("###MESSAGES###",messages)

    with open(outputfilename, 'w') as myfile:
        myfile.write(template)


def writeHTMLOutputPairs(loop_list, discrepancy_list, match_types, annotations, filename):
    """
    Write the list of loops in an HTML format that also shows
    the coordinate window and a heat map of all-against-all distances.
    Include information about matches to Matlab motif groups.
    num is the group number
    """

    JS1 = '  <script src="./js/JSmol.min.nojq.js"></script>'
    JS2 = '  <script src="./js/jquery.jmolTools.2loop.js"></script>'   # advance by 2 with Next
    JS3 = '  <script src="./js/imagehandlinglocal.js"></script>'
    JS4 = '<script src="./js/jmolplugin.js" type="text/javascript"></script>'
    JS5 = '<script type="text/javascript" src="./js/heatmap.js"></script>'
    TEMPLATEPATH = './search/'
    OUTPUTPATH = "clusters/"

    print("Writing HTML file for %3d loops" % (len(loop_list)))


    pagetitle = "Pairs of loops"

    # build the table that lists the loops
    candidatelist = '<table style="white-space:nowrap;">\n'

    # write header line
    candidatelist += "<tr><th>S.</th><th>Show</th>"
    candidatelist += "<th>Loop ID</th>"
    candidatelist += "<th>Discrepancy</th>"
    candidatelist += "<th>Match Type</th>"
    candidatelist += "<th>Annotation</th>"
    candidatelist += "</tr>\n"

    # write one row for each candidate
    for i in range(0,len(loop_list)):
        pair_number = (i/2) + 1
        candidate = loop_list[i]
        discrepancy = discrepancy_list[i]
        match_type = match_types[i]
        annotation = annotations[i]

        candidatelist += '<tr><td>%d.</td><td><label><input type="checkbox" id="%d" class="jmolInline" data-coord="%s">&nbsp</label></td>' % (pair_number,i,candidate)
        candidatelist += "<td>%s</td>" % (candidate)

        if i % 2 == 0:
            candidatelist += "<td>%8.4f</td>" % (discrepancy)
            candidatelist += "<td>%s</td>" % (match_type)
        else:
            candidatelist += "<td></td>"
            candidatelist += "<td></td>"

        candidatelist += "<td>%s</td>" % (annotation)
        candidatelist += '</tr>\n'

    candidatelist += '</table>\n'


    # read template.html into one string
    with open(TEMPLATEPATH + 'template.html', 'r') as myfile:
        template = myfile.read()

    template = template.replace("###PAGETITLE###",pagetitle)
    template = template.replace("###sequencecolumn###","")
    template = template.replace("###QUERYNAME###","")
    template = template.replace("###SEEMODIFYQUERY###","")
    template = template.replace("###seeCSVOutput###","")
    template = template.replace("###CANDIDATELIST###",candidatelist)
    template = template.replace("###JS1###",JS1)
    template = template.replace("###JS2###",JS2)
    template = template.replace("###JS3###",JS3)
    template = template.replace("###JS4###",JS4)

    template = template.replace("###REFRESH###","")

    if False:
        template = template.replace("###JS5###",JS5)    # include heatmap.js code
        discrepancydata = '<script type="text/javascript">\n' + discrepancydata + '\n</script>'
        template = template.replace("###DISCREPANCYDATA###",discrepancydata)
    else:
        template = template.replace("###DISCREPANCYDATA###","")
        template = template.replace("###JS5###","")    # do not display a heat map


    template = template.replace("###MESSAGES###","")

    with open(filename, 'w') as myfile:
        myfile.write(template)

"""
writeHTMLOutputPairs(['IL_4V9F_001','IL_4V9F_002','IL_4V9F_003','IL_4V9F_004'],
    [0.343,0.343,0.2,0.2],["homologous", "", "geometric", ""],
    ["", "secondAnnotation", "thirdAnnoation", "final"], 'output/test_pairs.html')
"""
def order_group_by_discrepancy(loop_ids, loop_ids_to_discrepancy, centroid_loop):
    ordered_by_discrepancy_loop_ids = {}
    for loop_id in loop_ids:
        ordered_by_discrepancy_loop_ids[loop_id] = loop_ids_to_discrepancy[centroid_loop][loop_id]

    ordered_by_discrepancy_loop_ids = {k: v for k, v in sorted(ordered_by_discrepancy_loop_ids.items(), key=lambda item: item[1])}
    return ordered_by_discrepancy_loop_ids


def get_matrix_for_consensus_interactions(loop_id_to_core_units):
    # pairTypes = ['glycosidicBondOrientation','chiDegree','pairsStacks','BPh','BR','sO','crossingNumber']
    pairTypes = ['basepair','pairsStacks','BPh','BR','others'] # change names stacks
    all_bptypes = ['ncHW', 'cSS', 'cWW', 'cWH', 'cHH', 'cHS', 'ntHW', 'ntHH', 'cSH', 'ntWS', 'ncSS', 'ncWH', 'ntWH', 'cHW', 'cWS', 'tHS', 'ntSH', 'tSH', 'tSW', 'ncHH', 'ncSH', 'tHW', 'ntSW', 'tHH', 'cSW', 'ncHS', 'ncWW', 'ntWW', 'tWS', 'ntSS', 'ntHS', 'tWW', 'ncSW', 'tSS', 'tWH', 'ncWS']
    all_stacks= ['s35','s55','s33','s53','ns35', 'ns55', 'ns33']

    all_unit_pairs = {}
    all_unit_pairs['basepair'] = set()
    all_unit_pairs['pairsStacks'] = set()
    all_unit_pairs['BR'] = set()
    all_unit_pairs['BPh'] = set()
    all_unit_pairs['others'] = set()

    loop_id_type_pair_to_interaction = defaultdict()
    for loop_id, nts in loop_id_to_core_units.items():
        pdb_id = loop_id.split("_")[1]

        # next line must be slow; interactions are already loaded for search spaces
        fr3d_pair_to_interaction_list, fr3d_pair_to_crossing_number = get_fr3d_pair_to_interaction_list(pdb_id, near=True, standard=True)

        loop_id_type_pair_to_interaction[loop_id] = {pairType: {} for pairType in pairTypes}

        # some loops are closed with AU or CG but not even annotated as near; need to be cWW
        corePositionPair = '1 - ' + str(len(nts))
        loop_id_type_pair_to_interaction[loop_id]['basepair'][corePositionPair] = 'cWW'

        for oneUnit in range(len(nts)):
            for anotherUnit in range(oneUnit+1, len(nts)):
                pair = (nts[oneUnit], nts[anotherUnit])
                consensus_interaction = fr3d_pair_to_interaction_list[pair]
                if len(consensus_interaction) >= 1:
                    # name for the pair of core positions
                    corePositionPair = str(oneUnit + 1) + ' - ' + str(anotherUnit + 1)
                    for i in consensus_interaction:
                        if i == 'cp':
                            # coplanar, not something we are tracking
                            pass
                        elif "R" in i:
                            # cSR or tSR, not something we are tracking
                            pass
                        elif i in all_bptypes:
                            # regularize cWWa, cWw, and other things like that
                            i = i.replace('a','').replace("w","W").replace("h","H").replace("s","S")
                            loop_id_type_pair_to_interaction[loop_id]['basepair'][corePositionPair] = i
                            all_unit_pairs['basepair'].add(corePositionPair)
                        elif i.startswith('c') or i.startswith('t') or i.startswith('nc') or i.startswith('nt'):
                            # more robust than using a list
                            i = i.replace('a','').replace("w","W").replace("h","H").replace("s","S")
                            loop_id_type_pair_to_interaction[loop_id]['basepair'][corePositionPair] = i
                            all_unit_pairs['basepair'].add(corePositionPair)
                        elif i in all_stacks:
                            loop_id_type_pair_to_interaction[loop_id]['pairsStacks'][corePositionPair] = i
                            all_unit_pairs['pairsStacks'].add(corePositionPair)
                        elif 'BR' in i:
                            loop_id_type_pair_to_interaction[loop_id]['BR'][corePositionPair] = i
                            all_unit_pairs['BR'].add(corePositionPair)
                        elif 'BPh' in i:
                            loop_id_type_pair_to_interaction[loop_id]['BPh'][corePositionPair] = i
                            all_unit_pairs['BPh'].add(corePositionPair)
                        else:
                            loop_id_type_pair_to_interaction[loop_id]['others'][corePositionPair] = i
                            all_unit_pairs['others'].add(corePositionPair)

    # fill in places with no interaction with empty string
    for pairType in pairTypes:
        for loop_id, nts in loop_id_to_core_units.items():
            for corePositionPair in all_unit_pairs[pairType]:
                if corePositionPair not in loop_id_type_pair_to_interaction[loop_id][pairType]:
                    loop_id_type_pair_to_interaction[loop_id][pairType][corePositionPair] = ''

    # print(loop_id_type_pair_to_interaction)

    # sort by corePositionPairs
    for loop_id in loop_id_type_pair_to_interaction:
        for pairtype in loop_id_type_pair_to_interaction[loop_id]:
            loop_id_type_pair_to_interaction[loop_id][pairtype] = {key: loop_id_type_pair_to_interaction[loop_id][pairtype][key] for key in sorted(loop_id_type_pair_to_interaction[loop_id][pairtype])}

    return loop_id_type_pair_to_interaction


def find_bp_signature(queries, loop_id_to_core_units, centroid_id, consensus_interactions):

    consensus = consensus_interactions['basepair']

    strand_info = queries[centroid_id]["loop_info"]["strand"]
    id_and_its_strand = []
    for unit_id in loop_id_to_core_units[centroid_id]:
        for strand in strand_info:
            if unit_id in strand:
                id_and_its_strand.append([unit_id, strand_info.index(strand)])

    dim_bp_matrix = len(loop_id_to_core_units[centroid_id])

    bp_matrix = [[0 for _ in range(dim_bp_matrix)] for _ in range(dim_bp_matrix)]

    paired_units = set()
    not_null_interactions = set()
    for unit_pair in consensus:
        if consensus[unit_pair] != '':
            x, y = [int(i) for i in unit_pair.split(' - ')]
            x, y = x-1, y-1
            bp_matrix[x][y] = consensus[unit_pair]
            paired_units.add(x)
            paired_units.add(y)
            not_null_interactions.add(unit_pair)

    unpaired_units = set()
    for i in range(dim_bp_matrix):
        if i not in paired_units:
            unpaired_units.add(i)

    bpSignature = []

    try:
        # set the outer basepair, which should always be cWW
        first_cww = '1 - %s' % (dim_bp_matrix)
        if consensus[first_cww] == 'cWW' or consensus[first_cww] == 'ncWW':
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

    # Keep R and L for IL, but use F for "fixed" in other loop types
    if not centroid_id.startswith('IL_'):
        for i in bpSignature:
            if i == "L" or i == "R":
                bpSignature[bpSignature.index(i)] = "F"
    bpSignature = "-".join(bpSignature)

    return bpSignature


def writeCSVOutput(group_id, loop_ids, loop_ids_to_discrepancy, loop_id_to_core_units, centroid_loop, loop_id_type_pair_to_interaction, queries, output_dir, consensus_interactions):
    with open(output_dir + '/' + 'MotifList.csv', 'at') as file:
        writer = csv.writer(file)
        for loop_id in loop_ids:
            writer.writerow([loop_id, 'Group_%03d' % group_id])
    with open(output_dir + '/' + 'MutualDiscrepancy.csv', 'at') as file:
        writer = csv.writer(file)
        for loop_i_id, inner_dict in loop_ids_to_discrepancy.items():
            for loop_j_id, discrepancy in inner_dict.items():
                # if loop_i_id != loop_j_id:
                writer.writerow([loop_i_id, round(discrepancy, 4), loop_j_id])
    with open(output_dir + '/' + 'MotifPositions.csv', 'at') as file:
        writer = csv.writer(file)
        for loop_id, nts in loop_id_to_core_units.items():
            for nt in nts:
                writer.writerow(['Group_%03d' % group_id, loop_id, nt, nts.index(nt)+1])
    # MotifLoopOrder: group_id, loop id, discrepancy order, similarity order
    # reordered_loop_ids is from order_group_by_similarity
    # need centroid loop id to get discrepancy order

    ordered_by_discrepancy_loop_ids = order_group_by_discrepancy(loop_ids, loop_ids_to_discrepancy, centroid_loop)
    with open(output_dir + '/' + 'MotifLoopOrder.csv', 'at') as file:
        writer = csv.writer(file)
        for loop_id in ordered_by_discrepancy_loop_ids:
            writer.writerow(['Group_%03d' % group_id, loop_id, list(ordered_by_discrepancy_loop_ids).index(loop_id) + 1, loop_ids.index(loop_id) + 1])

    bpSignature = find_bp_signature(queries, loop_id_to_core_units, centroid_loop, consensus_interactions)
    with open(output_dir + '/' + 'MotifBpSignatures.csv', 'at') as file:
        writer = csv.writer(file)
        writer.writerow(['Group_%03d' % group_id, bpSignature])



