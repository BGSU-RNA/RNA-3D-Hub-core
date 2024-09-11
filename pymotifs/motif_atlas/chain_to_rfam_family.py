import os, csv, gzip, time, sys #, wget
sys.path.append('search')
from fr3d_configuration import DATAPATHATLAS
# need to get rid of wget
if sys.version_info[0] < 3:
    from urllib import urlretrieve as urlretrieve
else:
    from urllib.request import urlretrieve as urlretrieve


def open_and_format_file(file_str, url, delim):
    '''
    This function checks to see if the string specified file from rfam
    exists, if not then download the file from rfam. Once the file is
    present this program reads and returns the tab separated file from rfam
    '''
    file_content = None
    if not os.path.exists(file_str):
        gz = None
        if url == "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/":
            if file_str == "Rfam.pdb":
                # gz = wget.download("{}{}".format(url,"Rfam.pdb.gz")) #need to instantiate url
                # ungzipped_file = gzip.open("Rfam.pdb.gz", "r") # added
                gz = urlretrieve("{}{}{}".format(url, file_str, ".gz"), os.path.join(DATAPATHATLAS, "Rfam.pdb.gz"))
                ungzipped_file = gzip.open(os.path.join(DATAPATHATLAS, "Rfam.pdb.gz"), "r")
            elif file_str == "family.txt":
                # gz = wget.download("{}{}".format(url,"database_files/{}.gz".format(file_str)))
                # ungzipped_file = gzip.open("family.text.gz", "r") # added
                gz = urlretrieve("{}{}{}".format(url, "database_files/", file_str, ".gz"), os.path.join(DATAPATHATLAS, ))
                ungzipped_file = gzip.open(os.path.join(DATAPATHATLAS, "family.text.gz"), "r")
            # ungzipped_file = gzip.open("{}.gz".format(file_str), "rb") # original
            # ungzipped_file = gzip.open("family.text.gz", "r")
            file_content = csv.reader(ungzipped_file.readlines(), delimiter=delim)
        elif url == "http://rna.bgsu.edu/rna3dhub/nrlist/download/3.248/all/csv":
            # equivalence_class_file = wget.download("{}".format(url))
            # equivalence_class_csv = open(equivalence_class_file, "r")
            equivalence_class_file = urlretrieve("{}".format(url), os.path.join(DATAPATHATLAS,"nr_list.3.248.csv"))
            equivalence_class_csv = open(os.path.join(DATAPATHATLAS, "nr_list_3.248.csv"), "r")
            file_content = csv.reader(equivalence_class_csv, delimiter=delim) # for each row in equivalence_class_csv:: index 0: format equivalence class id , index 1: represeantative IFE, index 2: all members
    else:
        text_file = open(file_str, "r")
        file_content = csv.reader(text_file.readlines(), delimiter=delim)
    return file_content
###


def map_chains_to_rfam_data(dict, rfam_pdb_list, rfam_family_list): # 4 spaces if this is blank
    '''
    This function maps chains to rfam data. Must run map_chains_to_equivalence_class
    for correct output. See example below.
    input: 1FFK|1|9 -> NR_all_25303.2, 4V9F|1|9
    output: 1FFK|1|9 -> NR_all_25303.2, 4V9F|1|9, RF0001, 5S ribosomal RNA, 1, 121
    '''
    # 1st. Looop over the pdb list and add the rfam family
    # and start and end nodes to the dictionary
    chain_to_rfamily_id = {}  # chain -> RF0001, chain start, chain end
    for line in rfam_pdb_list[1:]:
        chain = "{}|1|{}".format(line[1].upper(), line[2])  # ignore model number
        rfam_family = line[0]
        while chain in chain_to_rfamily_id.keys():
            chain += "*"
        chain_to_rfamily_id[chain] = rfam_family, line[3], line[4]
    # 2nd. Loop over the family list
    rfam_family_to_family = {}  # RF0001 -> 5S ribosomal RNA
    for line in rfam_family_list:
        family = line[-1]
        rfam_family = line[0]
        rfam_family_to_family[rfam_family] = family
    # 3rd. combine chain_to_rfamily_id and rfam_family_to_family by rfamily id
    rfamilyid_to_rfam_data = {}  # chain -> RF0001, chain start, chain end, 5S ribosomal RNA
    for key in chain_to_rfamily_id.keys():
        rfam_family = chain_to_rfamily_id[key][0]
        pdb_start = chain_to_rfamily_id[key][1]
        pdb_end = chain_to_rfamily_id[key][2]
        family = rfam_family_to_family[rfam_family]
        rfamilyid_to_rfam_data[key] = "{}\t{}\t{}\t{}\t".format(rfam_family, family, pdb_start, pdb_end)
    # 4th. join the rfam data to our dictionary dict
    for key in rfamilyid_to_rfam_data.keys():
        if key in dict.keys():
            dict[key] += rfamilyid_to_rfam_data[key]
        elif key.split("*")[0] in dict.keys():
            dict[key] = "{}".format(rfamilyid_to_rfam_data[key]) #add chain that maps to more than one rfam family as its own key
        else: #no equivalence class data in this case so add 2 tabs for alignment
            dict[key] = "{}".format(rfamilyid_to_rfam_data[key])
    for key in dict:
        if "+" in key:
            ifes = key.split("+")
            ife_in_ec_and_rfam = True
            for ife in ifes:
                if ife not in rfamilyid_to_rfam_data.keys():
                    ife_in_ec_and_rfam = False
            if ife_in_ec_and_rfam:
                rfamily_id = ""
                family_name = ""
                for ife in ifes:
                    rfamily_id += "{} + ".format(rfamilyid_to_rfam_data[ife].split("\t")[0])
                    family_name += "{} + ".format(rfamilyid_to_rfam_data[ife].split("\t")[1])
                rfamily_id = rfamily_id[:-2]  # remove plus sign
                family_name = family_name[:-2]
                dict[key] += "{}\t{}\t".format(rfamily_id, family_name)

# this used to be global, Adam put it into a function
def update_chain_to_rfamily_csv():

        # Urls
        rfam_url = "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/"  # current rfam build url

        # opening files
        try:
            rfam_family_txt = open_and_format_file("family.txt", rfam_url, "\t")
            rfam_pdb = open_and_format_file("Rfam.pdb", rfam_url, "\t")
        except:
            # print("something went wrong updating chain_to_fram_family.csv")
            # print("please manually update it")
            raise UserWarning("something went wrong updating chain_to_fram_family.csv\nplease manually update it")

        # reading files into individual arrays
        rfam_family_txt_arr = [line[:4] for line in rfam_family_txt]
        rfam_pdb_arr = [line for line in rfam_pdb]

        # data structure that we will use to map chains to their rfam and compound information.
        annotated_dictionary = {}

        # Next we add rfam data to this dictionary
        map_chains_to_rfam_data(annotated_dictionary, rfam_pdb_arr, rfam_family_txt_arr)
        # print(annotated_dictionary) # check the output of the function here
        # print(len(annotated_dictionary.keys())) # check to make sure all keys are present (there will likely be duplicates when there are chains mapped to multiple rfam families)

        with open('chain_to_rfam_family.csv', 'w') as csvfile:
            for key in annotated_dictionary.keys():
                csvfile.write("{}\t{}\n".format(key, annotated_dictionary[key]))
            print("csv write complete")


#########################################################################################
#### Past this point used to be a separate script called "proof of concept rfam interactions"
#########################################################################################

def read_equiv_class_csv_into_dict(newest_version = "current", break_ifes = False, cutoff = "all"):
    '''
    returns a dictionary of dict[chain]['equivalence class'] = 'NR_all_29909.1'
                            dict[chain]['quality rank'] = quality_rank = 1
    '''

    month_and_year = time.strftime("%b_%Y", time.gmtime()) # ex: "Sep_2023"

    # if not os.path.isfile("chain_to_rfam_family%s.csv".format(month_and_year)):
    #     update_chain_to_rfamily_csv()

    member_to_equiv_class = {}

    # by adding month and year suffix to filename, it will update monthly
    # filename = "nrlist_%s_all.csv" % (month_and_year)
    filename = "nrlist_%s.csv" % cutoff
    path_and_filename = os.path.join(DATAPATHATLAS, filename)

    # if not os.path.isfile(path_and_filename):
    if not False:
        print("Downloading " + filename + " to " + path_and_filename)
        try:
            urlretrieve("http://rna.bgsu.edu/rna3dhub/nrlist/download/" + str(newest_version) + "/all/csv", filename = path_and_filename)
        except:
            print("Download of equivalence class csv failed. Check internet connection.")
            return({})



    with open(path_and_filename) as file:
        rows = csv.reader(file, delimiter = ",")

        # for line in file:
        for row in rows:
            equiv_class = row[0]
            representative = row[1] # not used
            members = row[2].split(",")
            if break_ifes == False:
                for index, member in enumerate(members):
                    member_to_equiv_class[member] = {}
                    member_to_equiv_class[member]['equivalence class'] = equiv_class
                    member_to_equiv_class[member]['quality rank'] = index + 1
            if break_ifes == True:
                for index, member in enumerate(members):
                    for chain in member.split("+"):
                        member_to_equiv_class[chain] = {}
                        member_to_equiv_class[chain]['equivalence class'] = equiv_class
                        member_to_equiv_class[chain]['quality rank'] = index + 1


    return member_to_equiv_class



# member_to_equiv_class = read_equiv_class_csv_into_dict(newest_version = 3.257)

def read_chain_to_rfam():
    chain_to_rfam = {}

    with open("chain_to_rfamily.csv") as file:
        reading = csv.reader(file, delimiter="\t")
        for row in reading:
            # dict[chain] = [rfam ID, common name, nt start number, nt end number]
            chain_to_rfam[row[0]] = row[1:5]
    return chain_to_rfam


# for key, value in chain_to_rfam.items():
#   print(len(loops_and_strands))
    # print("{}: {}".format(key, value))

def which_rfam(loop):
    '''
    finding rfam family for each strand of this loop
    '''
    # chain_to_rfam[chain] example of return:
    # ['RF00001', '5S ribosomal RNA', '1', '119']

    filtered_chain_to_rfam = {}

    for strand in loop['strand']:
        # pulls the chain out of the first unit_id of each strand
        chain = "|".join(strand[0].split("|", 3)[0:3])

        try: # if chain in chain_to_rfam:
            family, molecule, start, stop = chain_to_rfam[chain]
            filtered_chain_to_rfam[chain] = {"family": family, "molecule": molecule}

        except KeyError: #else:
            filtered_chain_to_rfam[chain] = {"family": "no rfam family"}

    return(filtered_chain_to_rfam)


def which_rfam_and_equiv(loop, chain_to_rfam, member_to_equiv_class):
    '''
    finding rfam family and equivalence class for each strand of this loop
    '''
    # chain_to_rfam[chain] example of return:
    # ['RF00001', '5S ribosomal RNA', '1', '119']

    filtered_chain_to_rfam = {}

    for strand in loop['strand']:
        # pulls the chain out of the first unit_id of each strand
        chain = "|".join(strand[0].split("|", 3)[0:3])

        # make base keys
        try: # if chain in chain_to_rfam:
            family, molecule, start, stop = chain_to_rfam[chain]
            filtered_chain_to_rfam[chain] = {"family": family, "molecule": molecule}

        except KeyError: #else:
            filtered_chain_to_rfam[chain] = {"family": "no rfam family"}

        # add new equivalence keys
        try: # if chain in member_to_equiv_class:
            filtered_chain_to_rfam[chain]["equiv_class"] = member_to_equiv_class[chain]

        except KeyError: #else
            filtered_chain_to_rfam[chain]["equiv_class"] = "no equivalence class"

    return(filtered_chain_to_rfam)

def do_rfam_pools_match(loop1, loop2):
    # loops are expected to have the structure
    # loop[chain] = {"family": 'RF00001', "molecule": '5S ribosomal RNA'}

    rfam_pool1 = []

    rfam_pool2 = []

    for chain in loop1.values():
        rfam_pool1.append(chain['family'])

    for chain in loop2.values():
        rfam_pool2.append(chain['family'])

    if(set(rfam_pool1) == set(rfam_pool2)):
        return(True)
    else:
        return(False)


def do_equiv_pools_match(loop1, loop2):
    # loops are expected to have the structure
    # loop[chain] = {"family": 'RF00001', "molecule": '5S ribosomal RNA'}

    equiv_pool1 = []

    equiv_pool2 = []

    for chain in loop1.values():
        equiv_pool1.append(chain["equiv_class"])

    for chain in loop2.values():
        equiv_pool2.append(chain["equiv_class"])

    if(set(equiv_pool1) == set(equiv_pool2)):
        return(True)
    else:
        return(False)


def rfam_pool(loop):
    rfam_pool = []

    for chain in loop.values():
        rfam_pool.append(chain['family'])

    return(rfam_pool)


def equiv_pool(loop):
    rfam_pool = []

    for chain in loop.values():
        rfam_pool.append(chain['equiv_class'])

    return(rfam_pool)


def organize_equiv_per_chain(loop):
    '''
    this function is just to reorganize some parts
    of a dictionary for writing to csv
    '''
    output_dict = {}

    for chain in loop:
        output_dict[chain] = loop[chain]['family']

    return(output_dict)


def organize_family_per_chain(loop):
    '''
    this function is just to reorganize some parts
    of a dictionary for writing to csv
    '''
    output_dict = {}

    for chain in loop:
        output_dict[chain] = loop[chain]['equiv_class']

    return(output_dict)