"""

About

"""

__author__ = 'Anton Petrov'

import pdb, csv, logging, re, sys
import urllib2
from aMLSqlAlchemyClasses import session, PdbInfo


class PdbInfoLoader():
    """We need to check all files weekly in case anything changed"""

    def __init__(self):
        pdbs = []

    def get_all_rna_pdbs(self):
        """Gets a list of all rna-containing pdb files, including hybrids."""
        logging.info('Getting a list of all rna-containing pbds')
        url  = 'http://www.rcsb.org/pdb/rest/search'
        queryText = """
        <orgPdbQuery>
        <queryType>org.pdb.query.simple.ChainTypeQuery</queryType>
        <containsProtein>I</containsProtein>
        <containsDna>I</containsDna>
        <containsRna>Y</containsRna>
        <containsHybrid>I</containsHybrid>
        </orgPdbQuery>
        """
        req = urllib2.Request(url, data=queryText)
        f = urllib2.urlopen(req)
        result = f.read()

        if result:
            self.pdbs = result.split('\n')
            """filter out possible garbage in query results"""
            self.pdbs = filter(lambda x: len(x) == 4, self.pdbs)
            logging.info("Found number of PDB entries: %i", len(self.pdbs))
        else:
            logging.critical("Failed to retrieve results")
            sys.exit(2)

    def create_report_for_matlab(self):
        """
        """
#         for row in session.query(PdbInfo.structureId, PdbInfo.structureTitle, \
#                                  PdbInfo.experimentalTechnique, PdbInfo.releaseData, \
#                                  PdbInfo.structureAuthor, PdbInfo.resolution,
#                                  PdbInfo.source, func.group
        pass
# SELECT structureId,`structureTitle`,`experimentalTechnique`,`releaseDate`,`structureAuthor`, `resolution`,`source`,`chainLength`,
# group_concat(source ORDER BY chainLength DESC) AS source_
# FROM `pdb_info_copy`
# WHERE `entityMacromoleculeType` LIKE '%RNA%'
# AND structureId='2HGP'
# GROUP BY structureId
# keywords missing

    def _get_custom_report(self, pdb_id):
        """Gets a custom report in csv format for a single pdb file. Each chain
           is described in a separate line"""
        custom_report = '&customReportColumns=structureId,chainId,structureTitle,experimentalTechnique,depositionDate,releaseDate,revisionDate,ndbId,resolution,classification,structureMolecularWeight,macromoleculeType,structureAuthor,entityId,sequence,chainLength,db_id,db_name,molecularWeight,secondaryStructure,entityMacromoleculeType,ligandId,ligandIdImage,ligandMolecularWeight,ligandFormula,ligandName,ligandSmiles,InChI,InChIKey,hetId,Ki,Kd,EC50,IC50,deltaG,deltaH,deltaS,Ka,compound,plasmid,source,taxonomyId,biologicalProcess,cellularComponent,molecularFunction,ecNo,expressionHost,cathId,cathDescription,scopId,scopDomain,scopFold,pfamAccession,pfamId,pfamDescription,crystallizationMethod,crystallizationTempK,phValue,densityMatthews,densityPercentSol,pdbxDetails,unitCellAngleAlpha,unitCellAngleBeta,unitCellAngleGamma,spaceGroup,lengthOfUnitCellLatticeA,lengthOfUnitCellLatticeB,lengthOfUnitCellLatticeC,Z_PDB,rObserved,rAll,rWork,rFree,refinementResolution,highResolutionLimit,reflectionsForRefinement,structureDeterminationMethod,conformerId,selectionCriteria,fieldStrength,manufacturer,model,contents,solventSystem,ionicStrength,ph,pressure,pressureUnits,temperature,softwareAuthor,softwareName,version,method,details,conformerSelectionCriteria,totalConformersCalculated,totalConformersSubmitted,emdbId,emResolution,aggregationState,symmetryType,reconstructionMethod,specimenType&format=csv'
        url = 'http://www.rcsb.org/pdb/rest/customReport?pdbids=%s%s' % (pdb_id, custom_report)
        logging.info('Getting custom report for %s', pdb_id)
        f = urllib2.urlopen(url)
        result = f.read()
        if result:
            logging.info("Retrieved custom report for %s", pdb_id)
        else:
            logging.critical("Failed to retrieve results")
            sys.exit(2)
        lines = result.split('<br />');
        return lines

    def _compare_instances(self, A, B):
        """Compares two instances of the PdbInfo class. Reports any differences
           A - new instance, B - old instance"""
        for k,v in A.__dict__.iteritems():
            if k[0] == '_':
                continue # skip internal attributes
            if str(v) != str(getattr(B,k)): # compare as strings
                logging.warning('Structure %s, chain %s was updated', A.structureId, A.chainId)
                logging.warning('%s was: %s, became: %s', k, getattr(B,k), v)

    def _load_into_db(self, lines):
        """Compares the custom report data with what's already in the database,
           reports any discrepancies and stores the most recent version."""

        logging.info('Loading custom report')
        description = lines.pop(0) # contains field names
        keys = description.split(',')

        for line in lines:
            """one line per chain"""
            if len(line) < 10:
                continue # skip empty lines
            """replace quotechars that are not preceded and followed by commas,
            except for the beginning and the end of the string
            example: in 1HXL unescaped doublequotes in the details field"""
            line = re.sub('(?<!^)(?<!,)"(?!,)(?!$)',"'",line)
            """parse the line using csv reader"""
            reader = csv.reader([line],delimiter=',',quotechar='"')
            P = PdbInfo()
            for read in reader:
                for i,part in enumerate(read):
                    if part == '': # to save as NULL in the db
                        part = None
                    setattr(P,keys[i],part)

            logging.info('%s %s', P.structureId, P.chainId)
            """check if this pdb is present in the db"""
            existingP = session.query(PdbInfo). \
                                filter(PdbInfo.structureId==P.structureId). \
                                filter(PdbInfo.chainId==P.chainId).first()

            if existingP is not None and existingP != P:
                self._compare_instances(P, existingP)

            session.merge(P)
            session.commit()
        logging.info('Custom report saved in the database')

    def update_rna_containing_pdbs(self):
        """
        """
        try:
            if not pdbs:
                self.get_all_rna_pdbs()

            for pdb_id in self.pdbs:
                report = self._get_custom_report(pdb_id)
                self._load_into_db(report)
            logging.info('Successful update of RNA-containing pdbs')
            logging.info('%s', '+'*40)
        except:
            logging.critical('Update FAILED')
            sys.exit(2)



def usage():
    print __doc__


def main(argv):
    """
    """
    logging.basicConfig(level=logging.DEBUG)

    P = PdbInfoLoader()
    P.get_all_rna_pdbs()
#     P.pdbs = ['1HLX','1S72','2AVY']
    P.update_rna_containing_pdbs()


if __name__ == "__main__":
    main(sys.argv[1:])