"""
Query PDB for chain level information.

This asks PDB for the chain level information like sequence and compound for
the database. This works across several structures at once. This will merge the
data into the database if run several times the same data.
"""

import requests
from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils.pdb import CustomReportHelper

from pymotifs.pdbs.loader import Loader as PdbLoader


class Loader(core.SimpleLoader):
    merge_data = True
    dependencies = set([PdbLoader])
    @property
    def table(self):
        return mod.ChainInfo

    names = {
        'structureId': 'pdb_id',
        'chainId': 'chain_name',
        'classification': 'classification',
        'macromoleculeType': 'macromolecule_type',
        'entityId': 'entity_name',
        'sequence': 'sequence',
        'chainLength': 'chain_length',
        'source': 'source',
        'taxonomyId': 'taxonomy_id',
        'entityMacromoleculeType': 'entity_macromolecule_type',
        'compound': 'compound'
    }

    def rename(self, report):
        renamed = {}
        for key, name in self.names.items():
            renamed[name] = report.get(key)
        renamed['chain_id'] = report.get('chain_id')
        return renamed

    def get_ids(self, reports):
        pdbs = [report['structureId'] for report in reports]
        chains = [report['chainId'] for report in reports]

        mapping = {}
        known_pdbs = set()
        with self.session() as session:
            query = session.query(mod.ChainInfo.chain_id,
                                  mod.ChainInfo.pdb_id,
                                  mod.ChainInfo.chain_name).\
                filter(mod.ChainInfo.pdb_id.in_(pdbs)).\
                filter(mod.ChainInfo.chain_name.in_(chains))

            for result in query:
                mapping[(result.pdb_id, result.chain_name)] = result.chain_id
                known_pdbs.add(result.pdb_id)

        for report in reports:
            db_id = mapping.get((report['structureId'], report['chainId']))
            if not db_id and report['structureId'] in known_pdbs:
                self.logger.error("It seems a new chain %s was added to %s",
                                  report['chainId'], report['structureId'])

            if db_id:
                report['chain_id'] = db_id

        return reports

    def has_data(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(mod.ChainInfo).\
                filter_by(pdb_id=pdb).\
                limit(1)

            return bool(query.count())

    def query(self, session, pdb):
        """Generate a query to find all entries in chain_info for the given
        PDB id.  Added in November 2020 so this can be a SimpleLoader.
        Attributes
        ----------
        session : Session
            The `Session` to use.
        pdb : str
            The PDB id to use.
        Returns
        -------
        query : Query
            Returns an SqlAlchemy query for all entires in chain_info with
            the given PDB id.
        """
        return session.query(mod.ChainInfo).\
            filter_by(pdb_id=pdb)

    def olddata(self, pdbs, **kwargs):
        """ before November 2020, use PDB REST service
        Retrieve the data for each PDB file
        """

        helper = CustomReportHelper(fields=self.names)
        reports = helper(pdbs)
        data = []
        for report in self.get_ids(reports):
            renamed = self.rename(report)
            data.append(renamed)

        known = set(pdbs)
        seen = set(entry['pdb_id'] for entry in data)

        if len(seen) != len(known):
            missing = known - seen
            self.logger.error("Could not get info on all pdbs: %s",
                              ','.join(missing))

        return data

    def data(self, pdbs, **kwargs):
        """ after November 2020, use PDB graphQL to get data
        about each chain in the PDB file
        """

        # The GraphQL query is defined as a multi-line string.
        query = """
        {
          entry(entry_id:"XXXX") {
            entry {
              id
            }
            struct_keywords {
              pdbx_keywords
            }
            polymer_entities {
              rcsb_polymer_entity_container_identifiers {
                entry_id
                auth_asym_ids
                entity_id
              }
              entity_poly {
                rcsb_entity_polymer_type
                rcsb_sample_sequence_length
                pdbx_seq_one_letter_code_can
                type
              }
              rcsb_polymer_entity {
                pdbx_description
              }
              rcsb_entity_source_organism {
                ncbi_scientific_name
                ncbi_taxonomy_id
              }
            }
          }
        }
        """

        if isinstance(pdbs, str):
            pdbs = [pdbs]

        data = []    # will be a list over each pdb in pdbs and each chain in pdb

        for pdb in pdbs:

          currentquery = query.replace("XXXX",pdb)

          # the following line worked on rnatest but not on production
          #response = requests.post('http://data.rcsb.org/graphql', json={'query': currentquery})

          # here is a different way to do it
          currenturl = 'http://data.rcsb.org/graphql?query=' + currentquery
          response = requests.get(currenturl)

          if response.status_code == 200:
              result = response.json()
          else:
              self.logger.error("Could not get chain info for %s" % pdb)
              raise core.StageFailed("Could not load chain info for all pdbs")

          for chain_data in result["data"]["entry"]["polymer_entities"]:

            # some polymer_entities have more than one auth_asym_ids, which have different coords, see 7C7A
            for chain_id in chain_data["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"]:

              if not "ybrid" in chain_data["entity_poly"]["type"]:
                if "U" in chain_data["entity_poly"]["pdbx_seq_one_letter_code_can"]:
                    if "deoxy" in chain_data["entity_poly"]["type"]:
                      self.logger.info('DNA chain has U in it, type could be wrong')
                elif "T" in chain_data["entity_poly"]["pdbx_seq_one_letter_code_can"]:
                    if "olyribo" in chain_data["entity_poly"]["type"]:
                      self.logger.info('RNA chain has T in it, type could be wrong')

              renamed = {}
              renamed["pdb_id"]             = pdb
              renamed["chain_name"]         = chain_id
              renamed["entity_name"]        = chain_data["rcsb_polymer_entity_container_identifiers"]["entity_id"]
              renamed["classification"]     = result["data"]["entry"]["struct_keywords"]["pdbx_keywords"]
              renamed["macromolecule_type"] = chain_data["entity_poly"]["rcsb_entity_polymer_type"]
              renamed["sequence"]           = chain_data["entity_poly"]["pdbx_seq_one_letter_code_can"]
              renamed["chain_length"]       = chain_data["entity_poly"]["rcsb_sample_sequence_length"]
              renamed["entity_macromolecule_type"] = chain_data["entity_poly"]["type"]
              ######
              if chain_data["rcsb_entity_source_organism"] == None:
                chain_data["rcsb_entity_source_organism"] = [{'ncbi_scientific_name': None, 'ncbi_taxonomy_id': None}]
              renamed["taxonomy_id"]        = chain_data["rcsb_entity_source_organism"][0]["ncbi_taxonomy_id"]
              renamed["source"]             = chain_data["rcsb_entity_source_organism"][0]["ncbi_scientific_name"]              
              ######
              renamed["compound"]           = chain_data["rcsb_polymer_entity"]["pdbx_description"]

              data.append(renamed)

        return data