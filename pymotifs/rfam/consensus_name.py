"""
For each chain in each PDB file, try to assign a consensus name to the chain.
Do this by looking up mappings to Rfam families
Also by processing the pdb chain description.
Store information in the chain_property_value table
property can be:
    standardized_name
    source
    rfam_family
As of end of 2023, only store in chain_property_value when source, Rfam family, and standardized_name are known
"""

from cgi import print_arguments
import itertools as it
from operator import contains
import pymotifs.core as core
from pymotifs import models as mod
from pymotifs.chains.taxid_species_domain import Loader as DomainLoader
from pymotifs.rfam.map_to_rfam import Loader as MapLoader

class Loader(core.Loader):
    merge_data = True
    mark = False

    # needs to be done after the rfam mapping stage
    dependencies = set([DomainLoader,MapLoader])

    def to_process(self, pdbs, **kwargs):
        """
        Ignore list of pdb ids passed in.
        Find pdb ids that have an RNA chain in chain_info but that are not
        in chain_property_value, and pass those back.
        This way, we will work through the pdb ids that are not yet in
        chain_property_value.
        """

        # this could work better if it operated at the chain level
        # then if some chains in a PDB file are already mapped, it could still map more

        # pdb ids with an RNA chain in chain_info
        with self.session() as session:
            query = session.query(mod.ChainInfo.pdb_id).\
                filter(mod.ChainInfo.entity_macromolecule_type.contains('Polyribonucleotide'))
        not_in_table_pdbs = set([r.pdb_id for r in query])

        # pdb ids with a chain mapped to an Rfam family in chain_property_value
        with self.session() as session:
            query = session.query(mod.ChainPropertyValue.pdb_id).\
                filter(mod.ChainPropertyValue.property == 'rfam_family')
        in_table_pdbs = set([r.pdb_id for r in query])

        # return a list of the list of pdb ids that are not in chain_property_value
        return [sorted(not_in_table_pdbs - in_table_pdbs)]

    def has_data(self, pdb_dict, **kwargs):
        """
        This method can query the database after to_process but before data method
        to see if there is already data in the database, and if so, it returns True,
        and the data method does not have to work on this item.
        """

        return False

    def remove(self, pdb, **kwargs):
        """
        We need a method that would remove data from the database.
        But apparently it does not actually need to remove anything.
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


    def open_and_format_file(self,file_str, delim):
        """
        Generic function to read a csv file
        """
        import os, csv
        file_content = None
        text_file = open(file_str, "r")
        if os.path.exists(file_str):
            file_content = csv.reader(text_file.readlines(), delimiter=delim)
        return file_content


    # def tax_to_domain(self,tax_list):
    #     """
    #     Converts taxonomy flat file into a dictionary of tax ids to domain and lineage.
    #     """
    #     tax_dict = {}
    #     for row in tax_list:
    #         tax_id = row[0]
    #         domain = row[1]
    #         lineage = row[2]
    #         if tax_id not in tax_dict and domain != "":
    #             tax_dict[tax_id] = domain,lineage
    #     return tax_dict


    def get_taxid_to_domain(self):
        """
        Read mapping from taxonomy to domain into a dictionary
        """

        taxid_to_domain = {}
        with self.session() as session:
            query = session.query(mod.TaxidSpeciesDomain.taxonomy_id,
                                  mod.TaxidSpeciesDomain.domain)
            for row in query:
                if row.domain:
                    taxid_to_domain[row.taxonomy_id] = row.domain

        return taxid_to_domain


    def map_chains_to_rfam_data(self, rfam_pdb_list, rfam_family_list, pdb_id):
        """
        Create a dictionary mapping chains in pdb_id to rfam family information
        """

        # 1st. Loop over the pdb list and add the rfam family
        # and start and end nodes to the dictionary
        chain_to_rfamily_id = {}
        for line in rfam_pdb_list[1:]:
            if pdb_id.upper() in line[1].upper():
                chain = line[2]
                rfam_family = line[0]
                if not chain in chain_to_rfamily_id:
                    # family, bit_score, chain_start, chain_end
                    chain_to_rfamily_id[chain] = rfam_family, line[5], line[3], line[4]

        # 2nd. Loop over the family list
        # map Rfam id to Rfam family name
        # RF00001 -> 5S ribosomal RNA
        rfam_family_to_family = {}
        for line in rfam_family_list:
            family = line[-1]
            rfam_family = line[0]
            rfam_family_to_family[rfam_family] = family

        # 3rd. combine chain_to_rfamily_id and rfam_family_to_family by rfamily id
        # chain -> RF00001, chain start, chain end, 5S ribosomal RNA
        chain_to_rfam_data = {}
        for key in chain_to_rfamily_id:
            rfam_family = chain_to_rfamily_id[key][0]
            bit_score = chain_to_rfamily_id[key][1]
            chain_start = chain_to_rfamily_id[key][2]
            chain_end = chain_to_rfamily_id[key][3]
            family = rfam_family_to_family[rfam_family]
            chain_to_rfam_data[key] = [rfam_family, family, bit_score, chain_start, chain_end] # change dictionary to chain to rfam data
            #print ("%s | %s" % (key, chain_to_rfam_data[key]))

        return chain_to_rfam_data # keys are chains in the form pdbid|1|chain


    def domain_checking(self, rna_chains, all_chains, pdb_id):
        """
        Determine the biological source, which could be
        domain like archaea, bacteria, eukarya
        organelle like mitochondria, chloroplast
        other like virus and unknown

        Inputs:
        rna_chains is a list of dictionaries
        all_chains is protein chains and any other chains in the pdb file
        pdb_id
        """

        if len(rna_chains) == 0:
            return rna_chains

        checks = set()                        # collect pieces of evidence for a source
        from_title = False                    # note where the evidence comes from
        from_rfam_ncbi_mismatch = False
        from_ssu_lsu_ncbi_mismatch  = False
        clue_count_dict = {}

        voted_domain = None

        # check PDB title, which is stored with each chain
        title = rna_chains[0]['title']
        if "mitochond".upper() in title.upper() or "mitoriboso".upper() in title.upper():
            checks.add("Mitochondria")
            from_title = True
        elif "chloroplast".upper() in title.upper():
            checks.add("Chloroplast")
            from_title = True

        # check compound name for non-RNA chains since protein chains are often marked better
        mitochondria_count = 0.0  # use float so Python 2.7 does not do integer division
        chloroplast_count = 0.0
        indicator_count = 0.0
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

        # check if the domain of the organism is Eukarya but the chain maps to bacteria, mitochondria, or chlororplast ribosome
        ssu = None
        lsu = None
        lsu_rfam_id_set = {"RF02541","RF02546"}   # Rfam LSU families that mitochondria / chloroplast map into
        ssu_rfam_id_set = {"RF00177","RF02545"}   # Rfam SSU families that mitochondria / chloroplast map into
        for chain in rna_chains:
            if chain['tax_id'].isdigit() and 'rfam_family' in chain:
                if chain['domain'].startswith("Eukary"):
                    if chain['rfam_family'] in ssu_rfam_id_set:
                        ssu = chain
                    elif chain['rfam_family'] in lsu_rfam_id_set:
                        lsu = chain
        # look for both SSU and LSU mapping to an organelle somehow
        if ssu and lsu:
            from_ssu_lsu_ncbi_mismatch = True

        # check mismatch between rfam domain and ncbi domain   [This check is at the chain level]
        for chain in rna_chains:
            from_rfam_ncbi_mismatch = False
            lsu_ncbi_mismatch = False
            ssu_ncbi_mismatch = False
            from_ssu_lsu_ncbi_mismatch  = False
            if chain['tax_id'].isdigit() and 'rfam_family' in chain:
                if chain['domain'].startswith("Eukary"):
                    if chain['rfam_family'] in lsu_rfam_id_set:
                        lsu_ncbi_mismatch = True
                    elif chain['rfam_family'] in ssu_rfam_id_set:
                        ssu_ncbi_mismatch = True
            if lsu_ncbi_mismatch or ssu_ncbi_mismatch:
                from_rfam_ncbi_mismatch = True

            # Gather checks into a dictionary
            other_pdb_count = len(all_chains)
            if len(all_chains) < 1:
                self.logger.info("PDB %s does not have non-RNA chains" % (pdb_id))
                other_pdb_count = 1
            clue_checks_dict = {"PDB_Title" :from_title, "PDB_Name": voted_domain, "Other_chains": (indicator_count / other_pdb_count > 0), "NCBI_Rfam_mismatch":from_rfam_ncbi_mismatch, "SSU_LSU_mismatch":from_ssu_lsu_ncbi_mismatch}
            row_counter = 0
            for clue in clue_checks_dict.values():
                if clue:
                    row_counter += 1
            clue_count_dict[chain['chain_name']] = [row_counter, checks, clue_checks_dict]    #chains are keys

        # check the domain source of each row at the pdb id level
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
                    self.logger.info("Error Mitochondrial and Chloroplast found in same chain of %s" % pdb_id)

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

        return rna_chains


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
                    #     # check manual domain with the consensus domain
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


    def data(self, pdb_ids, **kwargs):
        """
        Input is an entire list of pdb ids, so we only have to read data files once.
        Add information about RNA chains in this pdb file.
        This method gets called on each item in the list returned by the
        to_process method, one after another.
        We will get one pdb identifier at a time.
        """

        # Load data files
        self.logger.info("Reading data files")

        # NCBI Taxonomy files
        # taxonomy_file = "/usr/local/pipeline/hub-core/aux/consensus_naming_aux_files/tax_dictionary.csv"
        # tax_dict = self.tax_to_domain(self.open_and_format_file(taxonomy_file,'\t')) # rename to tax_to_dict using find and replace

        # read database table with mappings of taxonomy ids to species and domain
        taxid_to_domain = self.get_taxid_to_domain()

        # Rfam files
        # rfam_pdb_file = "/usr/local/pipeline/hub-core/aux/consensus_naming_aux_files/Rfam.pdb"
        rfam_pdb_file = "/usr/local/pipeline/alignments/pdb_chain_to_rfam.txt"

        # load an Rfam data file listing all Rfam families and information about the family
        rfam_family_file = "/usr/local/pipeline/hub-core/aux/consensus_naming_aux_files/family.txt"
        rfam_family_txt = self.open_and_format_file(rfam_family_file, '\t') # returns a list of list of strings
        rfam_family_txt_arr = [line[:4] for line in rfam_family_txt]

        # read mappings of chains to Rfam families
        # TODO: need to be more sophisticated because more than one mapping may be listed for a chain
        #
        rfam_pdb = self.open_and_format_file(rfam_pdb_file, '\t')
        rfam_pdb_arr = [line for line in rfam_pdb if not line[0] == 'RF00000']

        # standard names file
        consensus_names_file = "/usr/local/pipeline/hub-core/aux/consensus_naming_aux_files/manual_consensus_names.tsv"
        consensus_names = self.open_and_format_file(consensus_names_file, '\t')
        consensus_arr = [line for line in consensus_names]
        rfam_family_to_standard_name = self.make_rfam_family_to_standard_name(consensus_arr)

        # store data about chains by their pdb id
        pdb_to_chains = {}

        # query the database for all chains in the chain_info table with these pdb_ids
        with self.session() as session:
            chain_query = session.query(mod.ChainInfo.pdb_id,
                                    mod.ChainInfo.chain_name,
                                    mod.ChainInfo.entity_macromolecule_type,
                                    mod.ChainInfo.chain_length,
                                    mod.ChainInfo.compound,
                                    mod.ChainInfo.taxonomy_id,
                                    mod.PdbInfo.title).\
                outerjoin(mod.PdbInfo, mod.ChainInfo.pdb_id == mod.PdbInfo.pdb_id).\
                filter(mod.ChainInfo.pdb_id.in_(pdb_ids))

            for r in chain_query: # r is a structured variable (like a row)

                if not r.pdb_id in pdb_to_chains:
                    pdb_to_chains[r.pdb_id] = {'rna_chains' : [], 'other_pdb_chains' : []}

                chain_dict = {}
                chain_dict['pdb_id'] = str(r.pdb_id)
                chain_dict['chain_name'] = str(r.chain_name)
                chain_dict['title'] = str(r.title)
                chain_dict['tax_id'] = str(r.taxonomy_id)
                chain_dict['compound'] = str(r.compound)
                chain_dict['chain_length'] = str(r.chain_length)
                # chain_dict['macromolecule'] = str(r.entity_macromolecule_type)

                # add domain according to taxonomy id
                chain_dict['domain'] = taxid_to_domain.get(chain_dict['tax_id'], "None")

                if 'polyribonucleotide' in r.entity_macromolecule_type.lower():
                    # this is an RNA chain
                    pdb_to_chains[r.pdb_id]['rna_chains'].append(chain_dict)
                else:
                    # not an RNA chain
                    pdb_to_chains[r.pdb_id]['other_pdb_chains'].append(chain_dict)

        data = [] # array containing annotated rna chains (domain, rfam family, standard name)

        # loop over pdb ids, process all chains for each one
        for pdb_id in sorted(pdb_to_chains.keys()):
            self.logger.info("Processing PDB file %s" %(pdb_id))

            # map chains to rfam data
            chain_to_rfam_data = self.map_chains_to_rfam_data(rfam_pdb_arr,rfam_family_txt_arr,pdb_id)

            rna_chains = pdb_to_chains[pdb_id]['rna_chains']
            other_pdb_chains = pdb_to_chains[pdb_id]['other_pdb_chains']

            #print(chain_to_rfam_data)

            if len(rna_chains) == 0:
                continue

            for chain_dict in rna_chains:
                # add rfam data
                if chain_dict['chain_name'] in chain_to_rfam_data:
                    chain_dict['rfam_family'] = chain_to_rfam_data[chain_dict['chain_name']][0]
                    chain_dict['rfam_description'] = chain_to_rfam_data[chain_dict['chain_name']][1]
                    chain_dict['rfam_bit_score'] = chain_to_rfam_data[chain_dict['chain_name']][2]

            # check for domain conflicts
            rna_chains = self.domain_checking(rna_chains, other_pdb_chains, pdb_id)

            # add standard name
            self.consensus_naming(rna_chains, rfam_family_to_standard_name)

            # Make dictionaries for this PDB id.  Each one becomes a row in chain_property_value table
            for chain in rna_chains:
                domain = chain.get('domain','')
                rfam_family = chain.get('rfam_family','')
                std_name = chain.get('std_name','')

                if domain == 'None':
                    domain = ''

                print("%s | %4s | %12s | %8s | %s" % (pdb_id,chain['chain_name'], domain, rfam_family, std_name))
                self.logger.info("%s | %4s | %12s | %8s | %s" % (pdb_id,chain['chain_name'], domain, rfam_family, std_name))

                # only save data if all fields are present; aim for completeness
                if domain and rfam_family and std_name:

                    # standard name dictionary
                    if 'std_name' in chain:
                        standard_name_dict = {}
                        standard_name_dict['pdb_id'] = pdb_id
                        standard_name_dict['chain'] = chain['chain_name']
                        standard_name_dict['property'] = "standardized_name"
                        standard_name_dict['value'] = chain['std_name']
                        data.append(standard_name_dict)

                    # Domain dictionary
                    if 'domain' in chain:
                        if not chain['domain'].startswith("None"):
                            domain_dict = {}
                            domain_dict['pdb_id'] = pdb_id
                            domain_dict['chain'] = chain['chain_name']
                            domain_dict['property'] = "source"
                            domain_dict['value'] = chain['domain']
                            data.append(domain_dict)

                    # Rfam family dictionary
                    if 'rfam_family' in chain:
                        rfam_family_dict = {}
                        rfam_family_dict['pdb_id'] = pdb_id
                        rfam_family_dict['chain'] = chain['chain_name']
                        rfam_family_dict['property'] = "rfam_family"
                        rfam_family_dict['value'] = chain['rfam_family']
                        data.append(rfam_family_dict)

        for row in data:
            #print(row)
            self.logger.info(row)

        self.logger.info('Ready to add %d pieces of data about chains' % len(data))
        print('Ready to add %d pieces of data about chains' % len(data))

        # temporarily do not write data to the database
        # data = []

        if len(data) > 0:
            # yield rows to the chain_property_value table
            for row in data:
                yield mod.ChainPropertyValue(**row)
        else:
            raise core.Skip("No chain name data found for pdb id %s" % (pdb_id))

