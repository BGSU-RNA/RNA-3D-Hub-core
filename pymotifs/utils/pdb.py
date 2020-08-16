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

    def __call__(self, dates=(None, None)):
        """Get a list of all rna-containing pdb files, including hybrids. Raise
        a specific error if it fails.

        :kwargs: Keyword arguments to define additional properties, such as
        date, to search by. See `.xml` for details.
        :returns: The response body if this succeeds.
        """

        try:
            fields = ['structureId', 'entityMacromoleculeType', 'releaseDate']
            helper = CustomReportHelper(fields=fields)
        except Exception as err:
            logger.exception(err)
            raise GetAllRnaPdbsError("Failed getting all PDBs")

        data = helper('*')
        ids = set()
        for entry in data:
            entity_type = entry['entityMacromoleculeType']
            if entity_type and 'RNA' in entity_type and \
                    self.within_date(dates, entry['releaseDate']):
                ids.add(entry['structureId'])
        return sorted(ids)


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
        logger.info('Found %d unique structures', len(unique))
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
