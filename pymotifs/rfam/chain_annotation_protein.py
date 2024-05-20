from pymotifs import core
from pymotifs import models as mod
import csv
import requests
import threading
import sys
import gzip
import time
from ftplib import FTP
import os


class Loader(core.SimpleLoader):

    def get_pdb_without_name(self, pdb):
        """ Query the database to find all PDB structures containing protein
        chains without name
        """

        with self.session() as session:
            query = session.query(mod.ChainInfo.pdb_id, mod.ChainInfo.chain_name) \
                           .filter(mod.ChainInfo.entity_macromolecule_type.contains('Polypeptide'))
        protein_pdbs = set([(r.pdb_id.lower(), r.chain_name) for r in query])

        with self.session() as session:
            query = session.query(mod.ChainPropertyValue)
        in_table_pdbs = set([(r.pdb_id.lower(), r.chain) for r in query])

        not_in_table = protein_pdbs - in_table_pdbs
        not_in_table_pdbs = set([entry[0] for entry in not_in_table])
        return list(not_in_table_pdbs)

    def get_protein_name(self, accession, unp_name_dict):
        uniprot_base_url = "https://rest.uniprot.org/uniprotkb/"
        uniprot_complete_url = uniprot_base_url + accession
        try:
            response = requests.get(uniprot_complete_url).json()
            if "recommendedName" in response["proteinDescription"]:
                protein_name = response["proteinDescription"]["recommendedName"]["fullName"]["value"]
                unp_name_dict[accession] = protein_name
        except Exception as e:
            print("Error getting protein name for accession {}: {}".format(accession, str(e)))

    def get_uniprot_naming(self, accession_list):
        unp_name_dict = {}
        threads = []
        for accession in accession_list:
            print("Getting name for the following accession", accession)
            t = threading.Thread(target=self.get_protein_name, args=(accession, unp_name_dict))
            threads.append(t)
            t.start()
        for t in threads:
            t.join()
        return unp_name_dict

    def get_pdb_unp_mapping(self, pdb_list):
        filepath = "/usr/local/pipeline/hub-core/aux/pdb_unp_mapping_files/pdb_chain_uniprot.csv.gz"
        unp_mapping_dict = {}
        unp_accessions = set()

        with gzip.open(filepath, 'rb') as f:
            content = cStringIO.StringIO(f.read())
            csv_reader = csv.reader(content, delimiter=",")
            next(csv_reader)
            next(csv_reader)
            for line in csv_reader:
                if line[0] in pdb_list:
                    identifier = (line[0], line[1])
                    print(identifier)
                    unp_mapping_dict[identifier] = line[2]
                    unp_accessions.add(line[2])

        return unp_mapping_dict, unp_accessions

    def download_pdb_mapping_file(self):
        # Define the FTP server credentials
        ftp_server = "ftp.ebi.ac.uk"

        # Define the remote file path on the FTP server
        remote_file_path = "/pub/databases/msd/sifts/csv/pdb_chain_uniprot.csv"

        # Define the local file path to save the downloaded file
        local_file_path = "/usr/local/pipeline/hub-core/aux/pdb_unp_mapping_files/pdb_chain_uniprot.csv"

        try:
            # Connect to the FTP server
            ftp = FTP(ftp_server)
            ftp.login()

            # Download the file from the FTP server
            with open(local_file_path, "wb") as local_file:
                ftp.retrbinary("RETR " + remote_file_path, local_file.write)

            # Compress the downloaded file using gzip
            with open(local_file_path, "rb") as f_in:
                with gzip.open(local_file_path + ".gz", "wb") as f_out:
                    f_out.writelines(f_in)

            # Delete the original uncompressed file
            os.remove(local_file_path)

            # Close the FTP connection
            ftp.quit()

            print("PDB-UniProt mapping file downloaded and compressed successfully!")

        except Exception as e:
            print("An error occurred while downloading and compressing the file:", str(e))

    def query(self, session, pdb):
        """The query method checks whether data is available in the db
        and if not, recomputes the data. I made the output of the query
        to be empty in order to recompute the data
        """

        chain_property_value = mod.ChainPropertyValue

        count = session.query(chain_property_value).\
                       filter(chain_property_value.pdb_id == None)

        return count

    def data(self, pdb, **kwargs):

        start_time = time.time()
        self.download_pdb_mapping_file()
        pdb_without_name = self.get_pdb_without_name(pdb)
        unp_mapping, unp_accessions = self.get_pdb_unp_mapping(pdb_without_name)
        unp_names = self.get_uniprot_naming(list(unp_accessions))

        print("Loading data into ChainPropertyValue table")
        for k, v in unp_mapping.iteritems():
            unp_name = unp_names.get(v, None)
            if unp_name:
                yield mod.ChainPropertyValue(pdb_id=k[0],
                                chain=k[1],
                                property="unp_accession",
                                value=v)
                yield mod.ChainPropertyValue(pdb_id=k[0],
                                chain=k[1],
                                property="unp_name",
                                value=unp_name)
        print("Done loading data into the table")

        end_time = time.time() - start_time
        print(end_time)