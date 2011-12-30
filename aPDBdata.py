import urllib2
import pdb, csv, logging, re, sys

from aMLSqlAlchemyClasses import *


class PdbInfoUpdater():
    """
    """

    def __init__(self):
        pdbs = []

    def get_all_rna_pdbs(self):
        """
        """
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
            logging.info("Found number of PDB entries: %i", result.count('\n'))
            self.pdbs = result.split('\n')
            """filter out possible garbage in query results"""
            self.pdbs = filter(lambda x: len(x) == 4, self.pdbs)
        else:
            logging.warning("Failed to retrieve results")
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

    def get_custom_report(self, pdb_id):
        """
        """
        custom_report = '&customReportColumns=structureId,chainId,structureTitle,experimentalTechnique,depositionDate,releaseDate,revisionDate,ndbId,resolution,classification,structureMolecularWeight,macromoleculeType,structureAuthor,entityId,sequence,chainLength,db_id,db_name,molecularWeight,secondaryStructure,entityMacromoleculeType,ligandId,ligandIdImage,ligandMolecularWeight,ligandFormula,ligandName,ligandSmiles,InChI,InChIKey,hetId,Ki,Kd,EC50,IC50,deltaG,deltaH,deltaS,Ka,compound,plasmid,source,taxonomyId,biologicalProcess,cellularComponent,molecularFunction,ecNo,expressionHost,cathId,cathDescription,scopId,scopDomain,scopFold,pfamAccession,pfamId,pfamDescription,crystallizationMethod,crystallizationTempK,phValue,densityMatthews,densityPercentSol,pdbxDetails,unitCellAngleAlpha,unitCellAngleBeta,unitCellAngleGamma,spaceGroup,lengthOfUnitCellLatticeA,lengthOfUnitCellLatticeB,lengthOfUnitCellLatticeC,Z_PDB,rObserved,rAll,rWork,rFree,refinementResolution,highResolutionLimit,reflectionsForRefinement,structureDeterminationMethod,conformerId,selectionCriteria,fieldStrength,manufacturer,model,contents,solventSystem,ionicStrength,ph,pressure,pressureUnits,temperature,softwareAuthor,softwareName,version,method,details,conformerSelectionCriteria,totalConformersCalculated,totalConformersSubmitted,emdbId,emResolution,aggregationState,symmetryType,reconstructionMethod,specimenType&format=csv'
        url = 'http://www.rcsb.org/pdb/rest/customReport?pdbids=%s%s' % (pdb_id, custom_report)
        f = urllib2.urlopen(url)
        result = f.read()
        if result:
            logging.info("Retrieved custom report for %s", pdb_id)
        else:
            logging.warning("Failed to retrieve results")
            sys.exit(2)

        lines = result.split('<br />');
        return lines


    def compare_instances(self, A, B):
        """
        """
        logging.warning('Structure %s, chain %s was updated', A.structureId, A.chainId)
        for k,v in A.__dict__.iteritems():
            if k[0] == '_':
                continue # skip internal attributes
            if str(v) != str(getattr(B,k)): # compare as strings
                logging.warning('%s was: %s, became: %s', k, v, getattr(B,k))


    def load_into_db(self, pdb_id):
        """
        """
        lines = self.get_custom_report(pdb_id)

        description = lines.pop(0)
        keys = description.split(',')

        for line in lines:
            if len(line) < 10:
                continue
            """replace quotechars that are not preceded and followed by commas,
            except for the beginning and the end of the string
            example: in 1HXL unescaped doublequotes in the details field"""
            line = re.sub('(?<!^)(?<!,)"(?!,)(?!$)',"'",line)
            reader = csv.reader([line],delimiter=',',quotechar='"')
            P = PdbInfo()
            for read in reader:
                for i,part in enumerate(read):
                    if part == '':
                        part = None
                    setattr(P,keys[i],part)

            logging.info('%s %s', P.structureId, P.chainId)
            existingP = session.query(PdbInfo). \
                                filter(PdbInfo.structureId==P.structureId). \
                                filter(PdbInfo.chainId==P.chainId).first()

            if existingP is not None and existingP != P:
                self.compare_instances(P, existingP)

            session.merge(P)
            session.commit()
#             pdb.set_trace()


    def update(self):
        """
        """
        self.get_all_rna_pdbs()
#         self.pdbs = ['1HLX']
        for pdb_id in self.pdbs:
            self.load_into_db(pdb_id)
        logging.info('SUCCESSFUL UPDATE')


logging.basicConfig(filename='pdbinfoupdate.log', filemode='w', level=logging.DEBUG)
P = PdbInfoUpdater()
try:
    P.update()
except:
    logging.warning('Update FAILED')
    sys.exit(2)
