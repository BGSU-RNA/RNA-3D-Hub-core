import re
import csv
import logging
import datetime

from pymotifs import utils

logger = logging.getLogger(__name__)


class GetAllRnaPdbsError(Exception):
    """Raised when when we cannot get a list of all pdbs."""
    pass


class RnaPdbsHelper(object):
    """A helper class to get a list of all RNA containing PDBS
    """

    url = 'http://www.rcsb.org/pdb/rest/search'

    def __init__(self):
        self.helper = utils.WebRequestHelper(method='post', parser=self.parse)

    def parse(self, response):
        return filter(lambda x: len(x) == 4, response.text.split("\n"))

    def __call__(self):
        """Get a list of all rna-containing pdb files, including hybrids. Raise
           a specific error if it fails.
        """

        logger.debug('Getting a list of all rna-containing pdbs')
        query_text = """
        <orgPdbQuery>
        <queryType>org.pdb.query.simple.ChainTypeQuery</queryType>
        <containsProtein>I</containsProtein>
        <containsDna>I</containsDna>
        <containsRna>Y</containsRna>
        <containsHybrid>I</containsHybrid>
        </orgPdbQuery>
        """
        headers = {
            'content-type': 'application/x-www-form-urlencoded'
        }

        try:
            return self.helper(self.url, headers=headers, data=query_text)
        except:
            logger.error("Could not get all rna containing pdbs")
            raise GetAllRnaPdbsError("Failed to get list of RNA PDBs")


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
        self.helper = utils.WebRequestHelper(method='get', parser=self.parse)
        if fields:
            self.fields = fields

    def parse(self, response):
        lines = response.text.split('<br />')
        description = lines.pop(0)
        keys = description.split(',')
        print(response.text)

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

    def __call__(self, pdb_id):
        """Gets a custom report in csv format for a single pdb file. Each chain
           is described in a separate line
        """

        params = {
            'customReportColumns': ','.join(self.fields),
            'format': 'csv',
            'pdbids': pdb_id
        }
        logger.info('Getting custom report for %s', pdb_id)
        return self.helper(self.url, params=params)


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
        return self.ftp.cwd('/pub/pdb/data/status/obsolete.dat',
                            parser=self.parse)
