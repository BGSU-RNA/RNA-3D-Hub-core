"""
Temporary stage to add modified nucleotides to loop_positions table

This module contains a loader to load all unit level information into the
database.
"""

from cgi import print_arguments
import itertools as it
from operator import contains
import pymotifs.core as core
from pymotifs import models as mod
from pymotifs.utils import units
from pymotifs import utils
# from pymotifs.download import Downloader
# from pymotifs.pdbs.info import Loader as PdbLoader
# from pymotifs.units import Loader as UnitInfoLoder
from sqlalchemy import and_
from collections import defaultdict
from Bio.Alphabet import ThreeLetterProtein
from sqlalchemy import desc
from sqlalchemy import asc
from sqlalchemy import not_




class Loader(core.Loader):
    merge_data = True
    mark = False

    # this stage is only run on existing chains, no dependencies needed
    dependencies = set([])

    def to_process(self, pdbs, **kwargs):
        '''
        After we get some files fixed and some not fixed, we may want to
        only work on the files that we didn't fix yet, by doing a query here
        and only returning the files that need work.
        For now, just return the list of pdbs that is passed into this function.
        '''
        with self.session() as session:
            # query = session.query(mod.ChainPropertyValue.pdb_id).\
            query = session.query(mod.ChainInfo.pdb_id).\
                filter(mod.ChainInfo.entity_macromolecule_type.contains('Polyribonucleotide'))
        not_in_table_pdbs = set([r.pdb_id for r in query])
        
        with self.session() as session:
            # query = session.query(mod.ChainPropertyValue.pdb_id).\
            query = session.query(mod.ChainPropertyValue.pdb_id)
        in_table_pdbs = set([r.pdb_id for r in query])
        # return sorted(set(pdbs) - in_table_pdbs)
        return(sorted(not_in_table_pdbs - in_table_pdbs))
        # return sorted(not_in_table_pdbs)

    def has_data(self, pdb_dict, **kwargs):
        """
        This method can query the database after to_process but before data method
        to see if there is already data in the database, and if so, it returns True,
        and the data method does not have to work on this item.
        """

        return False

    def remove(self, pdb, **kwargs):

        """
        with self.session() as session:
            query = session.query(mod.LoopInfo).filter_by(pdb_id=pdb)
            ids = [result.loop_id for result in query]

        if not ids:
            return True

        with self.session() as session:
            return session.query(mod.LoopPositions).\
                filter(mod.LoopPositions.loop_id.in_(ids)).\
                delete(synchronize_session=False)
        """

        return True

    # def current(self, unit_id): # current row info
    #     """Get the current data for the correspondence.
    #     """

    #     with self.session() as session:
    #         info = session.query(mod.UnitInfo).get(unit_id)
    #         return utils.row2dict(info)


    # def type_query(self, unit_ids, **kwargs):
    #     d = {}
    #     structure = self.structure(unit_ids[0][:4])
    #     # print(units.component_type(self.structure(unit_ids[:4]).unit_id()))
    #     # print(units.component_type(self.structure(unit_ids[:4]).residues('5Z9W|1|A|ALA|297||||P_35')))
    #     for base in structure.residues(polymeric=None):
    #         d.update({base.unit_id():units.component_type(base)})
    #     for unit_id in unit_ids:
    #         if unit_id.split('|')[3] in AA:
    #             d.update({unit_id:'aa'})
    #     return d

    '''
    This function checks to see if the string specified file from rfam 
    exists, if not then download the file from rfam. Once the file is 
    present this program reads and returns the tab separated file from rfam
    '''

    def open_and_format_file(self,file_str, delim):
        import os, csv
        file_content = None
        text_file = open(file_str, "r")
        if os.path.exists(file_str):
            file_content = csv.reader(text_file.readlines(), delimiter=delim)
        return file_content
    '''
    Converts taxonomy flat file into a dictionary of tax ids to domain and lineage.
    '''
    def tax_to_domain(self,tax_list):
        tax_dict = {}
        for row in tax_list:
            tax_id = row[0]
            domain = row[1]
            lineage = row[2]
            if tax_id not in tax_dict and domain != "":
                tax_dict[tax_id] = domain,lineage
        return tax_dict
    '''
    This function maps chains to rfam data.
    '''
    def map_chains_to_rfam_data(self,rfam_pdb_list, rfam_family_list, pdb_id): # Make this more efficient
        # 1st. Loop over the pdb list and add the rfam family
        # and start and end nodes to the dictionary
        # This function reads two data files (rfam.pdb, family.txt)
        chain_to_rfamily_id = {}  # chain -> RF0001, chain start, chain end
        for line in rfam_pdb_list[1:]:
            if pdb_id in line[1].upper():
                chain = line[2]  # ignore model number use the %s format 
                rfam_family = line[0]
                if chain not in chain_to_rfamily_id:
                    chain_to_rfamily_id[chain] = rfam_family, line[5], line[3], line[4]

        
        # 2nd. Loop over the family list
        rfam_family_to_family = {}  # RF0001 -> 5S ribosomal RNA
        for line in rfam_family_list:
            family = line[-1]
            rfam_family = line[0]
            rfam_family_to_family[rfam_family] = family
        
        # 3rd. combine chain_to_rfamily_id and rfam_family_to_family by rfamily id
        rfamilyid_to_rfam_data = {}  # chain -> RF0001, chain start, chain end, 5S ribosomal RNA
        for key in chain_to_rfamily_id:
            rfam_family = chain_to_rfamily_id[key][0]
            bit_score = chain_to_rfamily_id[key][1]
            chain_start = chain_to_rfamily_id[key][2]
            chain_end = chain_to_rfamily_id[key][3]
            family = rfam_family_to_family[rfam_family]
            rfamilyid_to_rfam_data[key] = [rfam_family, family, bit_score, chain_start, chain_end] # change dictionary to chain to rfam data
            # print ("%s | %s" % (key, rfamilyid_to_rfam_data[key]))

        return rfamilyid_to_rfam_data # keys are chains in the form pdbid|1|chain
    
    def domain_checking(self,rna_chains, all_chains, pdb_id):
        checks = set()
        from_title = False
        from_rfam_ncbi_mismatch = False
        from_ssu_lsu_ncbi_mismatch  = False
        clue_count_dict = {}
        
        mitochondria_count = 0.0
        chloroplast_count = 0.0
        indicator_count = 0.0
        voted_domain = None
        
        # checking title [This check is at the PDB level]
        # if len(rna_chains) > 0:
        title = rna_chains[0]['title']
        if "mitochond".upper() in title.upper() or "mitoriboso".upper() in title.upper():
            checks.add("Mitochondria")
            from_title = True
        elif "chloroplast".upper() in title.upper():
            checks.add("Chloroplast")
            from_title = True
        
        # checking compound name [This check is at the PDB level]
        for chain in all_chains:
            if "chloroplast".upper() in chain['compound'].upper():
                chloroplast_count += 1
                indicator_count += 1
            elif "mitochond".upper() in chain['compound'].upper() or "mitoriboso".upper() in chain['compound'].upper():
                mitochondria_count += 1
                indicator_count += 1

            m_percent =  100 * (mitochondria_count / len(all_chains))
            c_percent =  100 * (chloroplast_count / len(all_chains))
            threshold = 51
            if m_percent > threshold:
                voted_domain = "Mitochondria"
            elif c_percent > threshold:
                voted_domain = "Chloroplast"
        
        # checking if there is a mismatch between SSU chain and LSU chain domains [This check is at the PDB level]
        ssu = None
        lsu = None
        lsu_rfam_id_set = {"RF02541","RF02546"}
        ssu_rfam_id_set = {"RF00177","RF02545","RF02542"} 
        for row in rna_chains:
            if row['tax_id'].isdigit() and 'rfam_family' in row:
                if row['domain'].startswith("Eukarya"): # Should we check for mismatches outside of just Eukarya?
                    if row['rfam_family'] in ssu_rfam_id_set: 
                        ssu = row
                    elif row['rfam_family'] in lsu_rfam_id_set:
                        lsu = row
        if ssu and lsu:
            ssu_ncbi_domain = ssu['domain']
            lsu_ncbi_domain = lsu['domain']
            if ssu_ncbi_domain != lsu_ncbi_domain:
                from_ssu_lsu_ncbi_mismatch = True

        # checking mismatch between rfam domain and ncbi domain   [This check is at the chain level]

        for chain in rna_chains:
            from_rfam_ncbi_mismatch = False
            lsu_ncbi_mismatch = False
            ssu_ncbi_mismatch = False 
            from_ssu_lsu_ncbi_mismatch  = False
            if chain['tax_id'].isdigit() and 'rfam_family' in chain:
                if chain['domain'].startswith("Eukarya"):
                    if chain['rfam_family'] in lsu_rfam_id_set:
                        lsu_ncbi_mismatch = True
                    elif chain['rfam_family'] in ssu_rfam_id_set:
                        ssu_ncbi_mismatch = True
            if lsu_ncbi_mismatch or ssu_ncbi_mismatch:
                from_rfam_ncbi_mismatch = True

            #Gathering checks into a dictionary
            other_pdb_count = len(all_chains)
            if len(all_chains) < 1:
                self.logger.info("PDB %s is missing chains from all chains" % (pdb_id))
                other_pdb_count = 1
            clue_checks_dict = {"PDB_Title" :from_title, "PDB_Name": voted_domain, "Other_chains": (indicator_count / other_pdb_count > 0), "NCBI_Rfam_mismatch":from_rfam_ncbi_mismatch, "SSU_LSU_mismatch":from_ssu_lsu_ncbi_mismatch}
            row_counter = 0
            for clue in clue_checks_dict.values():
                if clue:
                    row_counter += 1
            clue_count_dict[chain['chain_name']] = [row_counter, checks, clue_checks_dict]    #chains are keys
        
        # checking the domain source of each row at the pdb id level
        for chain in rna_chains:
            source_str = "NCBI"
            # ith_row = rna_dict[pdb][i]
            consensus_domain = "None"
            if 'domain' in chain:
                consensus_domain = chain['domain'] #ith_row[7]
            num_clues = clue_count_dict[chain['chain_name']][0]
            check_set = clue_count_dict[chain['chain_name']][1]
            clue_booleans = clue_count_dict[chain['chain_name']][2]
            # Nicholas Reworking code stopped here February 15 @ 7:26PM
            if len(check_set) == 1:
                if num_clues > 1:
                    if clue_booleans["NCBI_Rfam_mismatch"]:
                        source_str = "Clue count (NCBI_Rfam_mismatch)"
                        consensus_domain = list(check_set)[0]
                    elif consensus_domain == "Eukarya":
                        source_str = "Clue count (NCBI = Eukarya and Mitochondrial / Chloroplast in PDB_Title/name/Other_chains)"
                        consensus_domain = list(check_set)[0]

                if chain['tax_id'].isdigit():
                    chain['source'] = source_str   
                else:
                    chain['source'] = "No_tax" 
                chain['clues'] = {'from_title':clue_booleans["PDB_Title"], 'from_compound_name':clue_booleans["PDB_Name"], 'percent_m_or_c_chains':clue_booleans["Other_chains"], 'from_Rfam_ncbi_domain_mismatch':clue_booleans["NCBI_Rfam_mismatch"], 'from_ssu_lsu_domain_mismatch':clue_booleans["SSU_LSU_mismatch"], 'num_clues': num_clues, 'check_set':check_set}
                chain['domain'] = consensus_domain

                if "Mitochondria" in checks and "Chloroplast" in checks:
                    print("Error Mitochondrial and Chloroplast found in same chain")
                    print(pdb_id)
                    print("")
        # check empty domain chains for clues as to their protien chains
        # for i in range(len(rna_dict[pdb])):
        domain_set = set()
        for chain in all_chains:
            if 'domain' in chain:
                chain_domain = chain['domain']
                domain_set.add(chain_domain)
        for chain in rna_chains:
            if 'domain' not in chain:
                if voted_domain is not None:
                    chain['domain'] = voted_domain
                    chain['source'] = "Majority of PDB chains have M or C in compound name"
                elif len(domain_set) == 1:
                        chain['domain'] = list(domain_set)[0]
                        chain['source'] = "Domain based on consensus PDB chains NCBI Taxonomy"

        # change the domain for each chain in a pdb id that has a mitochodrial or chloroplastic lsu or ssu chain
        mito_or_cholro_in_rna_chains = False
        new_source_str = ""
        for chain in rna_chains:
            if 'clues' in chain:
                if chain['clues']['from_Rfam_ncbi_domain_mismatch'] and chain['clues']['num_clues'] > 1: #17 is the mismatcj field
                    new_source_str = "PDB ID has Mitochondrial/chloroplastic LSU/SSU chain"
                    mito_or_cholro_in_rna_chains = True # if mito or chloro present in pdb id with an lsu or ssu chain and the # of clues > 1
        if mito_or_cholro_in_rna_chains:
            for chain in rna_chains:
                if len(chain['clues']['check_set']) == 1: #If there is a mito or chloro in chain
                    chain['domain'] = list(check_set)[0]
                    chain['source'] = new_source_str
                else:
                    chain['domain'] = voted_domain
                    chain['source'] = "Protein / RNA chains vote (compound name: Multiple domains found)"

        #return rna_dict.copy()
    
    def consensus_naming(self, rna_chains, rfam_family_to_standard_name):
        new_rfam = set()
        # for pdb in pdb_dict:
        for chain in rna_chains:
            pdb_compound_name = chain['compound']
            # rfam_id = chain['rfam_family']
            standard_name = None
            
            # chain['std_name'] = pdb_compound_name
            if 'rfam_family' in chain:
                if chain['rfam_family'] in rfam_family_to_standard_name :
                    chain['std_name'] = rfam_family_to_standard_name[chain['rfam_family']] # add consensus name
                    #     # checking manual domain with the consensus domain
                else:
                    # chain['std_name'] = None
                    if chain['rfam_family'] is not None:
                        if chain['rfam_family'].startswith("RF"):
                            new_rfam.add(chain['rfam_family'])
                # Archea rules          
                if "7S.S" in pdb_compound_name.upper(): # Rule: if the PDB compound contains "7S.S", say "Signal recognition particle domain S
                    standard_name = "Signal recognition particle domain S RNA; SRP S RNA"
                    chain['std_name'] = standard_name
                elif chain['tax_id'].isdigit():
                    if chain['domain'].startswith("Archaea") and '7S' in pdb_compound_name.upper(): # Rule:  if domain == 'Archaea' and '7S' in PDB_compound:
                        standard_name = "Signal recognition particle RNA; SRP RNA"
                        chain['std_name'] = standard_name
                    elif 'std_name' not in chain and chain['domain'].startswith("Archaea") and 'SRP' in pdb_compound_name.upper():
                        standard_name = "Signal recognition particle RNA; SRP RNA"
                        chain['std_name'] = standard_name
                    elif chain['domain'].startswith("Archaea") and chain['rfam_family'].startswith("RF00169"): # Rule when domain is archaea and it matches RF00169
                        standard_name = "Signal recognition particle RNA; SRP RNA"
                        chain['std_name'] = standard_name
                    elif 'std_name' not in chain and chain['domain'].startswith("Archaea") and '23S' in pdb_compound_name.upper() and int(chain['chain_length']) < 2000:
                        standard_name = "Large subunit ribosomal RNA fragment"
                        chain['std_name'] = standard_name              
                
                # RF00001 length rule: if chain length > 500 use PDB description as consensus name
                chain_length = chain['chain_length']
                if str(chain_length).isdigit():
                    if int(chain_length) > 500 and chain['rfam_family'] == "RF00001":
                        chain['std_name'] = pdb_compound_name
                if chain['tax_id'].isdigit():
                    if chain['domain'].startswith("Bacteria") and "16S" in pdb_compound_name.upper():
                        if int(chain_length) < 500:
                            standard_name = 'Small subunit ribosomal RNA fragment; SSU rRNA fragment'
                        else:
                            standard_name = 'Small subunit ribosomal RNA; SSU rRNA'
                        chain['std_name'] = standard_name   
                # if not chain['domain'] == 'Mitochondrial' and "16S" in pdb_compound_name and int(chain_length) < 500:
                #     name = 'Small subunit ribosomal RNA fragment; SSU rRNA fragment'

    def make_rfam_family_to_standard_name(self, consensus_name_list):
        rfam_family_to_standard_name = {}
        cn_list = consensus_name_list[1:]
        for line in cn_list:
            standard_name = line[3]
            rfam_family_to_standard_name[line[0]] = standard_name
        return rfam_family_to_standard_name

    def data(self, pdb_id, **kwargs):
        """
        This method gets called on each item in the list returned by the
        to_process method, one after another.
        We will get one pdb identifier at a time.
        """

        # query the database for all RNA chains in the Chain info table
        with self.session() as session:
            rna_query = session.query(mod.ChainInfo.pdb_id,
                                    mod.ChainInfo.chain_name,
                                    mod.ChainInfo.entity_macromolecule_type,
                                    mod.ChainInfo.chain_length,
                                    mod.ChainInfo.compound,
                                    mod.ChainInfo.taxonomy_id,
                                    mod.PdbInfo.title).\
                outerjoin(mod.PdbInfo, mod.ChainInfo.pdb_id == mod.PdbInfo.pdb_id).\
                filter(mod.ChainInfo.pdb_id == pdb_id, mod.ChainInfo.entity_macromolecule_type.contains('Polyribonucleotide')).\
                order_by(asc(mod.ChainInfo.chain_id))
        
        print("Current PDB: %s" %(pdb_id))
        data = [] # array containing annotated rna chains (domain, rfam family, standard name)

        if rna_query.count() > 0:
            
            print("Procesing PDB: %s" %(pdb_id))

            # query the database for data on all pdb chains that are not RNA in the Chain info table
            with self.session() as session:
                pdb_query = session.query(mod.ChainInfo.pdb_id,
                                        mod.ChainInfo.chain_name,
                                        mod.ChainInfo.entity_macromolecule_type,
                                        mod.ChainInfo.chain_length,
                                        mod.ChainInfo.compound,
                                        mod.ChainInfo.taxonomy_id,
                                        mod.PdbInfo.title).\
                    outerjoin(mod.PdbInfo, mod.ChainInfo.pdb_id == mod.PdbInfo.pdb_id).\
                    filter(mod.ChainInfo.pdb_id == pdb_id, not_(mod.ChainInfo.entity_macromolecule_type.contains('Polyribonucleotide'))).\
                    order_by(asc(mod.ChainInfo.chain_id))
            
            # Loading data files
            self.logger.info("Started Reading files!") # check timestamp in log file

            # NCBI Taxonomy files read
            taxonomy_file = "/usr/local/pipeline/hub-core/aux/consensus_naming_aux_files/tax_dictionary.csv"
            tax_dict = self.tax_to_domain(self.open_and_format_file(taxonomy_file,'\t')) # rename to tax_to_dict using find and replace
            
            # Rfam files read
            rfam_pdb_file = "/usr/local/pipeline/hub-core/aux/consensus_naming_aux_files/Rfam.pdb"
            rfam_family_file = "/usr/local/pipeline/hub-core/aux/consensus_naming_aux_files/family.txt"
            rfam_family_txt = self.open_and_format_file(rfam_family_file, '\t') # returns a list of list of strings
            rfam_pdb = self.open_and_format_file(rfam_pdb_file, '\t')
            rfam_family_txt_arr = [line[:4] for line in rfam_family_txt]
            rfam_pdb_arr = [line for line in rfam_pdb]
            rfam_dict = self.map_chains_to_rfam_data(rfam_pdb_arr,rfam_family_txt_arr,pdb_id)
            
            # standard names file read
            consensus_names_file = "/usr/local/pipeline/hub-core/aux/consensus_naming_aux_files/manual_consensus_names.tsv"
            consensus_names = self.open_and_format_file(consensus_names_file, '\t')
            consensus_arr = [line for line in consensus_names]
            self.logger.info("Ended Reading files!") # check time stamp in log file

            # Start processing PDB Files
            self.logger.info("Processing PDB file %s" %(pdb_id))

            other_pdb_chains = []
            rna_chains = []

            for r in rna_query: # r is a structured variable (like a row)
                chain_dict = {}
                chain_dict['pdb_id'] = str(r.pdb_id)
                chain_dict['chain_name'] = str(r.chain_name)
                # chain_dict['macromolecule'] = str(r.entity_macromolecule_type)
                chain_dict['chain_length'] = str(r.chain_length)
                chain_dict['compound'] = str(r.compound)
                chain_dict['title'] = str(r.title)
                chain_dict['tax_id'] = str(r.taxonomy_id)
                
                # adding ncbi taxonomy data
                if chain_dict['tax_id'] in tax_dict:
                    # chain_dict['lineage'] = tax_dict[chain_dict['tax_id']][1]
                    chain_dict['domain'] = tax_dict[chain_dict['tax_id']][0]
                else:
                    chain_dict['domain'] = "None"

                # adding rfam data
                if chain_dict['chain_name'] in rfam_dict:
                    chain_dict['rfam_family'] = rfam_dict[chain_dict['chain_name']][0]
                    chain_dict['rfam_description'] = rfam_dict[chain_dict['chain_name']][1]       
                    chain_dict['rfam_bit_score'] = rfam_dict[chain_dict['chain_name']][2]
                
                rna_chains.append(chain_dict)

            for r in pdb_query: # r is a structured variable (like a row)
                chain_dict = {}
                chain_dict['pdb_id'] = str(r.pdb_id)
                chain_dict['chain_name'] = str(r.chain_name)
                # chain_dict['macromolecule'] = str(r.entity_macromolecule_type)
                chain_dict['compound'] = str(r.compound)
                chain_dict['title'] = str(r.title)
                chain_dict['tax_id'] = str(r.taxonomy_id)
                if chain_dict['tax_id'] in tax_dict:
                    # chain_dict['lineage'] = tax_dict[chain_dict['tax_id']][1]
                    chain_dict['domain'] = tax_dict[chain_dict['tax_id']][0]
                other_pdb_chains.append(chain_dict)

            rfam_family_to_standard_name = self.make_rfam_family_to_standard_name(consensus_arr)


            # checking for domain conflicts
            self.domain_checking(rna_chains, other_pdb_chains, pdb_id) 

            # adding standard name
            self.consensus_naming(rna_chains, rfam_family_to_standard_name)

            rna_str = ""
            for line in rna_chains:
                if 'std_name' in line:
                    print("%s | %s | %s | %s\n" % (pdb_id,line['chain_name'], line['domain'], line['std_name']))
                    rna_str += "%s | %s | %s | %s\n" % (pdb_id, line['chain_name'], line['domain'], line['std_name'])

                else:
                    print("%s | %s | %s | %s\n" % (pdb_id,line['chain_name'], line['domain'], 'No standard name found'))                   
                    rna_str += "%s | %s | %s | %s\n" % (pdb_id, line['chain_name'], line['domain'], 'No standard name found')

            # write csv files here

            # rna_annotations = "/usr/local/pipeline/hub-core/aux/consensus_naming_aux_files/output.txt"
            # with open(rna_annotations, 'a') as f:
            #     f.write(rna_str)
            
            # Making dictionaries for PDB_ID
            for chain in rna_chains: #for chain_dict in other_pdb_chains. Check if its an rna chain, if so make this little dictionary if we have (name, doman rfam family), append it to the data list
                standard_name_dict = {}
                domain_dict = {}
                rfam_family_dict = {}
                chain_id = chain['chain_name']

                # standard name dict
                property_val = "standardized_name"
                if 'std_name' in chain:
                    value = chain['std_name']
                    standard_name_dict['pdb_id'] = pdb_id
                    standard_name_dict['chain'] = chain_id
                    standard_name_dict['property'] = property_val
                    standard_name_dict['value'] = value

                # Domain dict
                property_val = "source"
                if 'domain' in chain:
                    if not chain['domain'].startswith("None"):
                        value = chain['domain']
                        domain_dict['pdb_id'] = pdb_id
                        domain_dict['chain'] = chain_id
                        domain_dict['property'] = property_val
                        domain_dict['value'] = value
                
                # Rfam family dict
                property_val = "rfam_family"
                if 'rfam_family' in chain:
                    value = chain['rfam_family']
                    rfam_family_dict['pdb_id'] = pdb_id
                    rfam_family_dict['chain'] = chain_id
                    rfam_family_dict['property'] = property_val
                    rfam_family_dict['value'] = value

                # Appending dictionairies to the data list
                if len(standard_name_dict.values()) > 0:
                    data.append(standard_name_dict)
                    # print(standard_name_dict)
                    # self.logger.info(standard_name_dict)
                if len(domain_dict.values()) > 0:
                    # print(domain_dict)
                    # self.logger.info(domain_dict)
                    data.append(domain_dict)
                if len(rfam_family_dict.values()) > 0:
                    # print(rfam_family_dict)
                    # self.logger.info(rfam_family_dict)
                    data.append(rfam_family_dict)

        # yielding rows to the chain_property_value table if data is not empty
        if len(data) > 0:
            for row in data:
                # print("IT works!")
                yield mod.ChainPropertyValue(**row)
        else:
            raise core.Skip("RNA data not found for pdb id %s" % (pdb_id))
       
