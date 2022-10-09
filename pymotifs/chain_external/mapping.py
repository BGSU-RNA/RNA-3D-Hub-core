"""Stage to populate the units.info table.

This module contains a loader to load all unit level information into the
database.
"""

import itertools as it

import pymotifs.core as core

from pymotifs import models as mod
from pymotifs.utils import units
from pymotifs import utils
from sqlalchemy import and_
from collections import defaultdict
from Bio.Alphabet import ThreeLetterProtein
import requests
import os

AA = [seq.upper() for seq in ThreeLetterProtein().letters]

class Loader(core.Loader):
    merge_data = True
    mark = False
    allow_no_data = True


    dependencies = set([ChainInfo])
    """The dependencies for this stage."""


    def to_process(self, pdbs, **kwargs):
        # we used 5j7l for now
        pdbs = ['5j7l']
        return pdbs

    def remove(self, pdb_id, **kwargs):
        self.logger.info("Not removing anything")

    def has_data(self, unit_id, **kwargs):
        return 0

    def uniprot_chain_mapping(self,pdb,**kwargs):
        unique_uniprot_entries = set()
        chain_to_uniprot_info = {}
        uniprot_accession_info = {}
        base_url = "https://www.ebi.ac.uk/pdbe/graph-api/mappings/uniprot"
        uniprot_base_url = "https://rest.uniprot.org/uniprotkb/"

        # Get UniProt mappings for protein chains in a PDB structure
        complete_url = base_url + "/" + pdb
        # print('going to', complete_url)
        response = requests.get(complete_url).json()[pdb]['UniProt']
        # print(response)
        for accession in response:
            for entity in response[accession]['mappings']:
                chain = entity['chain_id']
                chain_to_uniprot_info[chain] = {"accession": accession, "organism": None, "recommended_name": None, "alternative_name": None}

        # create a set of UniProt accessions
        for chain, uniprot_info in chain_to_uniprot_info.items():
            unique_uniprot_entries.add(uniprot_info["accession"])

        # Get organism name, recommendedName and alternativeName for the UniProt accessions
        for accession in unique_uniprot_entries:
            uniprot_complete_url = uniprot_base_url + accession
            response = requests.get(uniprot_complete_url).json()
            organism_name = response["organism"]["scientificName"]
            recommended_protein_name = response["proteinDescription"]["recommendedName"]["fullName"]["value"]
            # there can be more than one alternative names but I'm just taking the first one
            alternative_protein_name = response["proteinDescription"]["alternativeNames"][0]["fullName"]["value"]
            uniprot_accession_info[accession] = {"organism": organism_name, "recommended_name": recommended_protein_name, "alternative_name": alternative_protein_name}


        for chain, uniprot_info in chain_to_uniprot_info.items():
            accession = uniprot_info["accession"]
            chain_to_uniprot_info[chain]["organism"] = uniprot_accession_info[accession]["organism"]
            chain_to_uniprot_info[chain]["recommended_name"] = uniprot_accession_info[accession]["recommended_name"]
            chain_to_uniprot_info[chain]["alternative_name"] = uniprot_accession_info[accession]["alternative_name"]

        return chain_to_uniprot_info

    def data(self, pdb, **kwargs):
        uniprot_chain_mapping_dict = self.uniprot_chain_mapping(pdb)
        for row in uniprot_chain_mapping_dict:
            chain_id = row
            # for info in uniprot_chain_mapping_dict[row]:
            external_value = uniprot_chain_mapping_dict[row]['accession']
            recommended_name = uniprot_chain_mapping_dict[row]['recommended_name']
            alternative_name = uniprot_chain_mapping_dict[row]['alternative_name']
            data = {
                    'chain_id': chain_id,
                    'pdb_id': pdb,
                    'external_value': external_value,
                    # 'recommended_name' = recommended_name
                    'external_property': alternative_name
                   }
            yield mod.UnitEcternalMapping(**data)
                





