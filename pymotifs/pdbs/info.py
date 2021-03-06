"""Update pdb_info table. This uses PDB RESTful services for getting custom
reports. This will fetch the data for all the given PDBs if any one does not
have data. The old data will be overwritten as needed.
"""

import requests
from pymotifs import core
from pymotifs.utils.pdb import CustomReportHelper

from pymotifs import models as mod


class Loader(core.SimpleLoader):
    """This loads information about entire files into our database.
    SimpleLoader does not insist on re-loading all PDB Info if one
    file is missing data, which is good after November 2020 when
    each file has to be queried separately from graphQL.
    """

    merge_data = True
    """We allow for merging data as the primary key is the pdb id."""

    dependencies = set()
    """This has no dependencies"""

    names = {
        'structureId': 'pdb_id',
        'structureTitle': 'title',
        'experimentalTechnique': 'experimental_technique',
        'depositionDate': 'deposition_date',
        'releaseDate': 'release_date',
        'revisionDate': 'revision_date',
        'ndbId': 'ndb_id',
        'resolution': 'resolution'
    }
    """This dict specifies the data to request from PDB's REST service (the
    keys) and what to call that data in our database (the values).
    """

    def rename(self, report):
        """Rename the columns of the report, and convert some values to the
        correct types. The columns will be renamed according to `Loader.names`
        and the resolution value will be converted to a float if it exists.

        Parameters
        ----------
        report : dict
            A dictonary with the keys in `self.names`.

        Returns
        -------
        renamed : dict
            The renamed dictonary
        """

        renamed = {}
        for key, name in self.names.items():
            renamed[name] = report.get(key)

        if renamed['resolution'] == '':
            renamed['resolution'] = None

        if renamed['resolution']:
            try:
                renamed['resolution'] = float(renamed['resolution'])
            except:
                renamed['resolution'] = None
                self.logger.error("Resoultion entry %s for %s is not a number",
                                  renamed['resolution'], renamed.get('pdb_id'))

        return renamed

    def has_data(self, pdb, **kwargs):
        """Check if we have data for the given PDB.

        Parameters
        ----------
        pdb : str
            The PDB to check

        Returns
        -------
        has_data : bool
            True if we have data for the PDB.
        """

        with self.session() as session:
            query = session.query(mod.PdbInfo).\
                filter_by(pdb_id=pdb).\
                limit(1)
            return bool(query.count())

    def query(self, session, pdb):
        """Generate a query to find all entries in PDBInfo for the given
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
            Returns an SqlAlchemy query for all entires in pdb_info with
            the given PDB id.
        """
        return session.query(mod.PdbInfo).\
            filter_by(pdb_id=pdb)

    def olddata(self, pdbs, **kwargs):
        """Compute a report of the given PDBs. This will use PDB's custom REST
        service to compute a report for the given PDBs. This will then be
        processed slightly for storage in our database. The columns used in the
        report will be the keys in `Loader.names`.

        Parameters
        ----------
        pdbs : list
            A list of PDB ids

        Returns
        -------
        reports : list
            A list of PdbInfo objects to write to the database.
        """

        helper = CustomReportHelper(fields=self.names.keys())
        data = helper(pdbs)

        print("pymotifs/pdbs/info.py gets this data:")
        print(data)
        print("for these pdb ids:")
        print(pdbs)

        if not data:
            raise core.StageFailed("Could not load data for all pdbs %s" %
                                   str(pdbs))

        if len(data) != len(pdbs):
            self.logger.error("Could not get all data for all pdbs")

        return [mod.PdbInfo(**self.rename(report)) for report in data]

    def data(self, pdbs, **kwargs):
        """New in November 2020.
        Get data from PDB's graphQL query.
        We loop over PDB IDs, make a query for each one, and
        run the query, then split out the data that is returned.

        Parameters
        ----------
        pdbs : list
            A list of PDB ids

        Returns
        -------
        reports : list
            A list of PdbInfo objects to write to the database.
        """

        # The GraphQL query is defined as a multi-line string.
        query = """
        {
          entry(entry_id:"XXXX") {
            entry {
              id
            }
            struct {
              title
            }
            exptl {
              method
            }
            rcsb_entry_info {
              resolution_combined
            }
            rcsb_accession_info {
              deposit_date
              initial_release_date
              revision_date
            }
            audit_author {
              name
            }
          }
        }
        """

        if isinstance(pdbs, str):
            pdbs = [pdbs]

        data = []

        for pdb in pdbs:

            self.logger.info("Using PDB graphQL to get data for %s file" % pdb)

            currentquery = query.replace("XXXX",pdb)

            # the following line worked on rnatest but not on production
            # response = requests.post('http://data.rcsb.org/graphql', json={'query': currentquery})

            # here is a different way to do it
            currenturl = 'http://data.rcsb.org/graphql?query=' + currentquery
            response = requests.get(currenturl)

            if response.status_code == 200:
                result = response.json()
            else:
                self.logger.error("Could not get PDB info for %s" % pdb)
                raise core.StageFailed("Could not load PDB info for all pdbs")

            renamed = {}
            renamed["pdb_id"] = pdb
            renamed["title"] = result["data"]["entry"]["struct"]["title"]
            renamed["experimental_technique"] = result["data"]["entry"]["exptl"][0]["method"]
            renamed["deposition_date"] = result["data"]["entry"]["rcsb_accession_info"]["deposit_date"][0:10]
            renamed["release_date"] = result["data"]["entry"]["rcsb_accession_info"]["initial_release_date"][0:10]
            renamed["revision_date"] = result["data"]["entry"]["rcsb_accession_info"]["revision_date"][0:10]
            renamed["ndb_id"] = pdb
            renamed["resolution"] = result["data"]["entry"]["rcsb_entry_info"]["resolution_combined"]
            renamed["authors"] = ", ".join([x["name"] for x in result["data"]["entry"]["audit_author"]])

            if renamed['resolution']:
                try:
                    renamed['resolution'] = float(renamed['resolution'][0])
                except:
                    renamed['resolution'] = None
                    self.logger.error("Resoultion entry for %s is not a number" % pdb)

            data.append(renamed)

        return [mod.PdbInfo(**report) for report in data]
