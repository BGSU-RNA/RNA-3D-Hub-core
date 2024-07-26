import re
import csv
import logging
import datetime
import requests

from pymotifs import utils

logger = logging.getLogger(__name__)

class GetAllRnaPdbsError(Exception):
    """Raised when when we cannot get a list of all pdbs."""
    pass


class RnaPdbsHelper(object):
    """
    A helper class to get a list of all RNA containing PDBs.
    Called by pymotifs/cli/setup.py
    """

    def within_date(self, dates, release):
        if not release:
            return None
        if not dates:
            return True

        if dates[0]:
            min_date = datetime.datetime.strptime(dates[0], r'%Y-%m-%d').date()
        else:
            min_date = datetime.date.min

        if dates[1]:
            max_date = datetime.datetime.strptime(dates[1], r'%Y-%m-%d').date()
        else:
            max_date = datetime.date.max

        release_date = datetime.datetime.strptime(release, r'%Y-%m-%d').date()

        return min_date <= release_date <= max_date

    def __call__(self, dates=(None, None), molecule_types = ["RNA","NA-hybrid"]):
        """
        Get a list of all rna-containing pdb files, including hybrids. Raise
        a specific error if it fails.

        :kwargs: Keyword arguments to define additional properties, such as
        date, to search by. See `.xml` for details.
        :returns: The response body if this succeeds.
        """

        # new code in November 2020 because the REST service is going away
        if dates[0]:
            earliest_date = dates[0]
        else:
            earliest_date = "1000-01-01"   # before PDB

        if dates[1]:
            latest_date = dates[1]
        else:
            latest_date = "3000-01-01"   # far in the future

        logger.info('Earliest date is %s, latest date is %s' % (earliest_date,latest_date))
        print('Earliest date is %s, latest date is %s' % (earliest_date,latest_date))

        resultIDs = []

        # change http to https, v1 to v2, pager to paginate on 2022-07-20
        # limit to 10,000 rows on 2023-08-09 when PDB imposed that limit
        # need to use the paginate option when the number of results exceeds 10,000
        url = "https://search.rcsb.org/rcsbsearch/v2/query?json=%7B%22query%22%3A%7B%22type%22%3A%22group%22%2C%22logical_operator%22%3A%22and%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22attribute%22%3A%22entity_poly.rcsb_entity_polymer_type%22%2C%22operator%22%3A%22exact_match%22%2C%22value%22%3A%22RNA%22%7D%7D%2C%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22operator%22%3A%22greater_or_equal%22%2C%22value%22%3A%221000-01-01T00%3A00%3A00Z%22%2C%22attribute%22%3A%22rcsb_accession_info.initial_release_date%22%7D%7D%2C%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22operator%22%3A%22less_or_equal%22%2C%22value%22%3A%223000-01-01T00%3A00%3A00Z%22%2C%22attribute%22%3A%22rcsb_accession_info.initial_release_date%22%7D%7D%5D%7D%2C%22request_options%22%3A%7B%22paginate%22%3A%7B%22start%22%3A0%2C%22rows%22%3A10000%7D%7D%2C%22return_type%22%3A%22entry%22%7D"
        url = url.replace("1000-01-01",earliest_date)
        url = url.replace("3000-01-01",latest_date)

        # query converted to regular text
        # never seems to work; maybe PDB simply won't accept these characters?
        # tried to convert this way, but doesn't work
        # import urllib.parse
        # print(urllib.parse.unquote(url))
        # url = 'https://search.rcsb.org/rcsbsearch/v2/query?json={"query":{"type":"group","logical_operator":"and","nodes":[{"type":"terminal","service":"text","parameters":{"attribute":"entity_poly.rcsb_entity_polymer_type","operator":"exact_match","value":"RNA"}},{"type":"terminal","service":"text","parameters":{"operator":"greater_or_equal","value":"1000-01-01T00:00:00Z","attribute":"rcsb_accession_info.initial_release_date"}},{"type":"terminal","service":"text","parameters":{"operator":"less_or_equal","value":"3000-01-01T00:00:00Z","attribute":"rcsb_accession_info.initial_release_date"}}]},"request_options":{"paginate":{"start":0,"rows":10000}},"return_type":"entry"}'
        # url = url.replace("1000-01-01",earliest_date)
        # url = url.replace("3000-01-01",latest_date)

        for molecule_type in molecule_types:
            currenturl = url.replace("RNA",molecule_type)
            # print("Querying PDB with url %s" % currenturl)
            logger.info("Querying PDB with url %s" % currenturl)
            response = requests.get(currenturl)
            jsonR = response.json()

            print("Number of %-9s rows received is %5d; limit is 10,000" % (molecule_type,len(jsonR["result_set"])))
            logger.info("Number of %-9s rows received is %5d; limit is 10,000" % (molecule_type,len(jsonR["result_set"])))

            for item in jsonR["result_set"]:
                resultIDs.append(item["identifier"])

        # except Exception as err:
        #     # for key in jsonR.keys():
        #     #     print(key,jsonR[key])
        #     print("PDB query error %s" % err)
        #     logger.exception(err)
        #     raise GetAllRnaPdbsError("Failed getting all PDBs")

        print("utils/pdb.py: Found %d distinct non-obsolete PDB ids of %s" % (len(set(resultIDs)),", ".join(molecule_types)))
        logger.info("Found %d distinct non-obsolete PDB ids of %s" % (len(set(resultIDs)),", ".join(molecule_types)))

        return(sorted(resultIDs))

        # code from before November 2020
"""
        try:
            fields = ['structureId', 'entityMacromoleculeType', 'releaseDate']
            helper = CustomReportHelper(fields=fields)
        except Exception as err:
            logger.exception(err)
            raise GetAllRnaPdbsError("Failed getting all PDBs")

        latestDate = '1000-01-01'

        data = helper('*')
        ids = set()
        for entry in data:
            entity_type = entry['entityMacromoleculeType']
            if entity_type and 'RNA' in entity_type and \
                    self.within_date(dates, entry['releaseDate']):
                ids.add(entry['structureId'])
                if entry['releaseDate'] > latestDate:
                	latestDate = entry['releaseDate']

        logger.info('Found %d unique PDB IDs', len(ids))
        logger.info('Most recent release date is %s', latestDate)

        print("IDs found by existing pipeline but not new search")
        print(sorted(ids.difference(set(resultIDs))))

        print("IDs found by new search but not existing pipeline")
        print(set(resultIDs).difference(ids))

        print("Old query found %d PDB ids" % len(ids))
        print("New query found %d PDB ids" % len(resultIDs))

        return sorted(ids)
"""

class CustomReportHelper(object):
    """A helper class to get a custom report from PDB.
    """

    url = 'http://www.rcsb.org/pdb/rest/customReport'

    fields = [
        'structureId', 'chainId', 'structureTitle',
        'experimentalTechnique', 'depositionDate', 'releaseDate',
        'revisionDate', 'ndbId', 'resolution', 'classification',
        'structureMolecularWeight', 'macromoleculeType',
        'structureAuthor', 'entityId', 'sequence', 'chainLength',
        'db_id', 'db_name', 'molecularWeight', 'secondaryStructure',
        'entityMacromoleculeType', 'hetId', 'Ki', 'Kd', 'EC50', 'IC50',
        'deltaG', 'deltaH', 'deltaS', 'Ka', 'compound', 'plasmid',
        'source', 'taxonomyId', 'biologicalProcess', 'cellularComponent',
        'molecularFunction', 'ecNo', 'expressionHost', 'cathId',
        'cathDescription', 'scopId', 'scopDomain', 'scopFold',
        'pfamAccession', 'pfamId', 'pfamDescription',
        'crystallizationMethod', 'crystallizationTempK', 'phValue',
        'densityMatthews', 'densityPercentSol', 'pdbxDetails',
        'unitCellAngleAlpha', 'unitCellAngleBeta', 'unitCellAngleGamma',
        'spaceGroup', 'lengthOfUnitCellLatticeA',
        'lengthOfUnitCellLatticeB', 'lengthOfUnitCellLatticeC', 'Z_PDB',
        'rObserved', 'rAll', 'rWork', 'rFree', 'refinementResolution',
        'highResolutionLimit', 'reflectionsForRefinement',
        'structureDeterminationMethod', 'conformerId', 'selectionCriteria',
        'contents', 'solventSystem', 'ionicStrength', 'ph', 'pressure',
        'pressureUnits', 'temperature', 'softwareAuthor', 'softwareName',
        'version', 'method', 'details', 'conformerSelectionCriteria',
        'totalConformersCalculated', 'totalConformersSubmitted', 'emdbId',
        'emResolution', 'aggregationState', 'symmetryType',
        'reconstructionMethod', 'specimenType'
    ]

    def __init__(self, fields=None):
        self.helper = utils.WebRequestHelper(method='post', parser=self.parse)
        if fields:
            self.fields = fields

    def parse(self, response):
        lines = response.text.split('<br />')
        description = lines.pop(0)
        keys = description.split(',')

        report = []
        for line in lines:
            """one line per chain"""
            if len(line) < 10:
                continue  # skip empty lines

            """replace quotechars that are not preceded and followed by commas,
            except for the beginning and the end of the string
            example: in 1HXL unescaped doublequotes in the details field"""

            line = re.sub('(?<!^)(?<!,)"(?!,)(?!$)', "'", line)
            reader = csv.reader([line], delimiter=',', quotechar='"')
            chain_dict = {}
            for read in reader:
                for i, part in enumerate(read):
                    if not part:
                        part = None  # to save as NULL in the db
                    chain_dict[keys[i]] = part
            report.append(chain_dict)
        return report

    def __unique__(self, reports):
        seen = set()
        unique = []
        for report in reports:
            entry = tuple(report.items())
            if entry not in seen:
                unique.append(report)
                seen.add(entry)
        return unique

    def __call__(self, pdb_id):
        """Gets a custom report in csv format for a single pdb file. Each chain
           is described in a separate line
        """
        ids = pdb_id
        if not isinstance(pdb_id, str):
            ids = ','.join(pdb_id)

        data = {
            'customReportColumns': ','.join(self.fields),
            'format': 'csv',
            'pdbids': ids
        }

        print("setup.py getting custom report")

        logger.info('Getting custom report for %s', pdb_id)
        result = self.helper(self.url, data=data)
        return self.__unique__(result)


class ObsoleteStructureHelper(object):
    def __init__(self):
        self.ftp = utils.FTPFetchHelper('ftp.wwpdb.org')

    def parse(self, text):
        data = []
        for line in text.split('\n'):
            if 'OBSLTE' in line:
                parts = line.split()
                data.append({
                    'id': parts[2],
                    'date': datetime.strptime(parts[1], '%d-%b-%y'),
                    'replaced_by': parts[3:]
                })
        return data

    def __call__(self):
        return self.ftp.action('/pub/pdb/data/status/obsolete.dat',
                               parser=self.parse)
