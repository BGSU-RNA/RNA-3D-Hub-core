import re
import csv
import logging
import datetime
import xml.etree.ElementTree as ET

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
        parser = lambda r: filter(lambda x: len(x) == 4, r.text.split("\n"))
        self.helper = utils.WebRequestHelper(method='post', parser=parser,
                                             allow_empty=True)

    def __date_refinement__(self, start, stop):
        refinement = {
            'queryType': 'org.pdb.query.simple.ReleaseDateQuery',
            'database_PDB_rev.date.comparator': 'between',
        }

        if start:
            refinement['database_PDB_rev.date.min'] = start.isoformat()

        if stop:
            refinement['database_PDB_rev.date.max'] = stop.isoformat()

        return refinement

    def __type_refinement__(self, rna='Y', dna='?', protein='?', hybrid='?'):
        return {
            'queryType': 'org.pdb.query.simple.ChainTypeQuery',
            'containsProtein': protein,
            'containsDna': dna,
            'containsRna': rna,
            'containsHybrid': hybrid
        }

    def xml(self, **kwargs):
        """Create the XML text to send to PDB for some search. This defaults to
        all RNA containing PDB files for any date. It does not limit if any
        other molecules are present. Given more than 1 of the other options it
        will and them all together.

        :contains: A dictonary containing any of rna, dna, protein and hybrid,
        where a value of True means the must contain that type of molecule,
        False meaning it should not and None (or not present) meaning ignore.
        Also if strings are given as values those will be used directly.
        :dates: A tuple of the start and end dates to limit the search to. If
        no start or end date is desired None, or some other falsey value should
        be used.
        :returns: A string of the xml data to send to PDB.
        """

        contains = {}
        for name, value in kwargs.get('contains', {}):
            if value is True:
                contains[name] = 'Y'
            elif value is False:
                contains[name] = 'N'
            elif value is None:
                contains[name] = '?'
            else:
                contains[name] = value

        refinements = [self.__type_refinement__(**contains)]
        if 'dates' in kwargs:
            start, stop = kwargs.get('dates')
            refinements.append(self.__date_refinement__(start, stop))

        root = ET.Element("orgPdbCompositeQuery", version='1.0')
        for index, refinement in enumerate(refinements):
            refine = ET.SubElement(root, "queryRefinement")
            level = ET.SubElement(refine, "queryRefinementLevel")
            level.text = str(index)
            if index != 0:
                ET.SubElement(refine, 'conjunctionType').text = 'and'

            refine_query = ET.SubElement(refine, "orgPdbQuery")
            for name, value in refinement.items():
                ET.SubElement(refine_query, name).text = value

        return ET.tostring(root)

    def __call__(self, **kwargs):
        """Get a list of all rna-containing pdb files, including hybrids. Raise
        a specific error if it fails.

        :kwargs: Keyword arguments to define additional properties, such as
        date, to search by. See `.xml` for details.
        :returns: The response body if this succeeds.
        """

        logger.debug('Getting a list of all rna-containing pdbs')
        query_text = self.xml(**kwargs)

        headers = {
            'content-type': 'application/x-www-form-urlencoded'
        }

        try:
            return self.helper(self.url, headers=headers, data=query_text)
        except Exception as err:
            logger.error("Could not get all rna containing pdbs")
            logger.exception(err)
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
        return self.ftp.cwd('/pub/pdb/data/status/obsolete.dat',
                            parser=self.parse)
