import requests
import time
import gzip
import csv
import io
from concurrent.futures import ThreadPoolExecutor
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from pymotifs import core
from pymotifs import models as mod
from sqlalchemy.exc import IntegrityError
import os


class Loader(core.SimpleLoader):

    PDB_MAPPING_FILE = "/usr/local/pipeline/hub-core/aux/pdb_unp_mapping_files/pdb_chain_uniprot.csv.gz"
    PDB_UNP_MAPPING_FILE_URL = "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/csv/pdb_chain_uniprot.csv"
    
    def configure_retry(self):
        """Set up retry logic for HTTP requests."""
        retries = Retry(
            total=5,  # Retry 5 times
            backoff_factor=1,  # Wait progressively longer between retries
            status_forcelist=[429, 500, 502, 503, 504],  # Retry for these HTTP errors
        )
        return HTTPAdapter(max_retries=retries)
    
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
        not_in_table_pdbs = set([(entry[0], entry[1]) for entry in not_in_table])
        return list(not_in_table_pdbs)
    
    def query(self, session, pdb):
        """The query method checks whether data is available in the db
        and if not, recomputes the data. I made the output of the query
        to be empty in order to recompute the data
        """

        chain_property_value = mod.ChainPropertyValue

        count = session.query(chain_property_value).\
                       filter(chain_property_value.pdb_id == None)

        return count

    def get_protein_name(self, accession, unp_name_dict):
        """Fetch protein name from UniProt for a given accession."""
        uniprot_base_url = "https://rest.uniprot.org/uniprotkb/"
        uniprot_complete_url = f"{uniprot_base_url}{accession}"

        session = requests.Session()
        session.mount("https://", self.configure_retry())

        try:
            response = session.get(uniprot_complete_url, timeout=10)
            response.raise_for_status()  # Raise an error for HTTP codes >= 400
            data = response.json()

            # Extract protein name
            if "recommendedName" in data.get("proteinDescription", {}):
                protein_name = data["proteinDescription"]["recommendedName"]["fullName"]["value"]
                unp_name_dict[accession] = protein_name
                print(f"Got protein name for accession {accession}")

            # Add a delay after a successful request
            time.sleep(0.5)

        except requests.exceptions.RequestException as e:
            print(f"Error getting protein name for accession {accession}: {e}")

    def get_uniprot_naming(self, accession_list):
        """Fetch UniProt names for a list of accessions with a thread pool."""
        unp_name_dict = {}

        def fetch_name(accession):
            self.get_protein_name(accession, unp_name_dict)

        with ThreadPoolExecutor(max_workers=10) as executor:  # Control parallelism
            executor.map(fetch_name, accession_list)

        return unp_name_dict

    def get_pdb_unp_mapping(self, pdb_list):
        """Parse the PDB-UniProt mapping file and extract mappings."""
        unp_mapping_dict = {}
        unp_accessions = set()

        # Convert pdb_list to a set for faster lookups
        pdb_set = set(pdb_list)

        with gzip.open(self.PDB_MAPPING_FILE, 'rb') as f:
            # Decode binary stream and read as CSV
            content = io.StringIO(f.read().decode('utf-8'))
            csv_reader = csv.reader(content, delimiter=",")

            # Skip two headers
            next(csv_reader)
            next(csv_reader)

            # Parse file and filter by PDB list tuples
            for line in csv_reader:
                identifier = (line[0], line[1])
                if identifier in pdb_set:
                    unp_mapping_dict[identifier] = line[2]
                    unp_accessions.add(line[2])

        return unp_mapping_dict, unp_accessions

    def download_pdb_mapping_file(self):
        """Download and compress the PDB-UniProt mapping file via HTTPS."""
        
        local_path = self.PDB_MAPPING_FILE.replace(".gz", "")

        try:
            # Download the file
            response = requests.get(self.PDB_UNP_MAPPING_FILE_URL, stream=True)
            response.raise_for_status()  # Raise an error for bad responses (4xx or 5xx)

            with open(local_path, "wb") as local_file:
                for chunk in response.iter_content(chunk_size=8192):
                    local_file.write(chunk)

            # Compress the file
            with open(local_path, "rb") as f_in, gzip.open(self.PDB_MAPPING_FILE, "wb") as f_out:
                f_out.writelines(f_in)

            os.remove(local_path)  # Remove the uncompressed file
            print("PDB-UniProt mapping file downloaded and compressed successfully!")

        except requests.exceptions.RequestException as e:
            print(f"Error downloading the file: {e}")

    def data(self, pdb, **kwargs):
        print("Downloading PDB mapping file...")
        self.download_pdb_mapping_file()
        print("Download is complete...")
        print("Parsing PDB-UniProt mappings...")
        pdb_without_name = self.get_pdb_without_name(pdb)
        print("How many PDBs have protein chains without name?", len(pdb_without_name))
        unp_mapping, unp_accessions = self.get_pdb_unp_mapping(pdb_without_name)

        print("How many UNP accesions?", len(unp_accessions))

        print("Fetching UniProt protein names...")
        unp_names = self.get_uniprot_naming(list(unp_accessions))

        print("Loading data into ChainPropertyValue table...")
        for (pdb_id, chain), unp_accession in unp_mapping.items():
            # Only process entries with valid UniProt names
            unp_name = unp_names.get(unp_accession)
            if not unp_name:
                print(f"Skipping entry for PDB: {pdb_id}, Chain: {chain}, UniProt: {unp_accession} (Name not available)")
                continue  # Skip entries without names

            print(f"Adding mapping information for PDB: {pdb_id}, Chain: {chain}, UniProt: {unp_accession}")
            try:
                yield mod.ChainPropertyValue(
                    pdb_id=pdb_id,
                    chain=chain,
                    property="unp_accession",
                    value=unp_accession,
                )
                yield mod.ChainPropertyValue(
                    pdb_id=pdb_id,
                    chain=chain,
                    property="unp_name",
                    value=unp_name,
                )
            except IntegrityError as e:
                # Log and skip duplicate entries
                print(f"Skipping duplicate entry for PDB: {pdb_id}, Chain: {chain}, Property: {e.params['property']}")
                continue

        print("Done loading data into the table.")