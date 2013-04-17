"""

Can be used in conjuction with nodetests. Just run
nodetests
from the project's home directory.

"""


import unittest

from PdbInfoLoader import PdbInfoLoader


class TestPdbInfoLoader(unittest.TestCase):

    ClassIsSetup = False

    def setUp(self):
        # If it was not setup yet, do it
        if not self.ClassIsSetup:
            print "Initializing testing environment"
            self.prepare()
            self.__class__.ClassIsSetup = True

    def prepare(self):
        unittest.TestCase.setUp(self)
        self.__class__.loader = PdbInfoLoader()
        self.__class__.loader.get_all_rna_pdbs()

    def test_length_of_rna_pdb_list(self):
        """Make sure that the pdb list is at least 2000 entries long"""
        self.assertTrue(len(self.loader.pdbs) > 2000)

    def test_pdb_list_contains_four_character_ids(self):
        """Make sure all ids are four characters long"""
        four_chars = filter(lambda x: len(x) == 4, self.loader.pdbs)
        self.assertTrue(len(self.loader.pdbs), len(four_chars))

    def test_get_custom_report(self):
        """Check that the report matches the standard"""
        line = self.loader._get_custom_report('1A60')
        result = ['structureId,chainId,structureTitle,experimentalTechnique,depositionDate,releaseDate,revisionDate,ndbId,resolution,classification,structureMolecularWeight,macromoleculeType,structureAuthor,entityId,sequence,chainLength,db_id,db_name,molecularWeight,secondaryStructure,entityMacromoleculeType,ligandId,ligandIdImage,ligandMolecularWeight,ligandFormula,ligandName,ligandSmiles,InChI,InChIKey,hetId,Ki,Kd,EC50,IC50,deltaG,deltaH,deltaS,Ka,compound,plasmid,source,taxonomyId,biologicalProcess,cellularComponent,molecularFunction,ecNo,expressionHost,cathId,cathDescription,scopId,scopDomain,scopFold,pfamAccession,pfamId,pfamDescription,crystallizationMethod,crystallizationTempK,phValue,densityMatthews,densityPercentSol,pdbxDetails,unitCellAngleAlpha,unitCellAngleBeta,unitCellAngleGamma,spaceGroup,lengthOfUnitCellLatticeA,lengthOfUnitCellLatticeB,lengthOfUnitCellLatticeC,Z_PDB,rObserved,rAll,rWork,rFree,refinementResolution,highResolutionLimit,reflectionsForRefinement,structureDeterminationMethod,conformerId,selectionCriteria,fieldStrength,manufacturer,model,contents,solventSystem,ionicStrength,ph,pressure,pressureUnits,temperature,softwareAuthor,softwareName,version,method,details,conformerSelectionCriteria,totalConformersCalculated,totalConformersSubmitted,emdbId,emResolution,aggregationState,symmetryType,reconstructionMethod,specimenType', '"1A60","A","NMR STRUCTURE OF A CLASSICAL PSEUDOKNOT: INTERPLAY OF SINGLE-AND DOUBLE-STRANDED RNA, 24 STRUCTURES","SOLUTION NMR","1998-03-04","1998-05-27","2003-04-01#2009-02-24","1A60","","RNA","13984.4","RNA","Kolk, M.H., Van Der Graaf, M., Wijmenga, S.S., Pleij, C.W.A., Heus, H.A., Hilbers, C.W.","1","GGGAGCUCAACUCUCCCCCCCUUUUCCGAGGGUCAUCGGAACCA","44","1A60","PDB","13984.4","GGGAGCUCAACUCUCCCCCCCUUUUCCGAGGGUCAUCGGAACCA#","Polyribonucleotide (RNA)","","","","","","","","","","","","","","","","","","TYMV PSEUDOKNOT","","Turnip yellow mosaic virus","12154","","","","","","","","","","","","","","","","","","","","90.0","90.0","90.0","P 1","1.0","1.0","1.0","1","","","","","","","","","","","750.0","BRUKER","AMX","","","","","","","","","X-PLOR","","TORSION ANGLE DYNAMICS/ SIMULATED ANNEALING","REFINEMENT DETAILS CAN BE FOUND IN THE JRNL CITATION ABOVE.","LEAST RESTRAINT VIOLATIONS","140","24","","","","","",""', '"1A60","A","NMR STRUCTURE OF A CLASSICAL PSEUDOKNOT: INTERPLAY OF SINGLE-AND DOUBLE-STRANDED RNA, 24 STRUCTURES","SOLUTION NMR","1998-03-04","1998-05-27","2003-04-01#2009-02-24","1A60","","RNA","13984.4","RNA","Kolk, M.H., Van Der Graaf, M., Wijmenga, S.S., Pleij, C.W.A., Heus, H.A., Hilbers, C.W.","1","GGGAGCUCAACUCUCCCCCCCUUUUCCGAGGGUCAUCGGAACCA","44","1A60","PDB","13984.4","GGGAGCUCAACUCUCCCCCCCUUUUCCGAGGGUCAUCGGAACCA#","Polyribonucleotide (RNA)","","","","","","","","","","","","","","","","","","TYMV PSEUDOKNOT","","Turnip yellow mosaic virus","12154","","","","","","","","","","","","","","","","","","","","90.0","90.0","90.0","P 1","1.0","1.0","1.0","1","","","","","","","","","","","600.0","BRUKER","AMX","","","","","","","","","X-PLOR","","TORSION ANGLE DYNAMICS/ SIMULATED ANNEALING","REFINEMENT DETAILS CAN BE FOUND IN THE JRNL CITATION ABOVE.","LEAST RESTRAINT VIOLATIONS","140","24","","","","","",""', '"1A60","A","NMR STRUCTURE OF A CLASSICAL PSEUDOKNOT: INTERPLAY OF SINGLE-AND DOUBLE-STRANDED RNA, 24 STRUCTURES","SOLUTION NMR","1998-03-04","1998-05-27","2003-04-01#2009-02-24","1A60","","RNA","13984.4","RNA","Kolk, M.H., Van Der Graaf, M., Wijmenga, S.S., Pleij, C.W.A., Heus, H.A., Hilbers, C.W.","1","GGGAGCUCAACUCUCCCCCCCUUUUCCGAGGGUCAUCGGAACCA","44","1A60","PDB","13984.4","GGGAGCUCAACUCUCCCCCCCUUUUCCGAGGGUCAUCGGAACCA#","Polyribonucleotide (RNA)","","","","","","","","","","","","","","","","","","TYMV PSEUDOKNOT","","Turnip yellow mosaic virus","12154","","","","","","","","","","","","","","","","","","","","90.0","90.0","90.0","P 1","1.0","1.0","1.0","1","","","","","","","","","","","500.0","BRUKER","AM","","","","","","","","","X-PLOR","","TORSION ANGLE DYNAMICS/ SIMULATED ANNEALING","REFINEMENT DETAILS CAN BE FOUND IN THE JRNL CITATION ABOVE.","LEAST RESTRAINT VIOLATIONS","140","24","","","","","",""', '"1A60","A","NMR STRUCTURE OF A CLASSICAL PSEUDOKNOT: INTERPLAY OF SINGLE-AND DOUBLE-STRANDED RNA, 24 STRUCTURES","SOLUTION NMR","1998-03-04","1998-05-27","2003-04-01#2009-02-24","1A60","","RNA","13984.4","RNA","Kolk, M.H., Van Der Graaf, M., Wijmenga, S.S., Pleij, C.W.A., Heus, H.A., Hilbers, C.W.","1","GGGAGCUCAACUCUCCCCCCCUUUUCCGAGGGUCAUCGGAACCA","44","1A60","PDB","13984.4","GGGAGCUCAACUCUCCCCCCCUUUUCCGAGGGUCAUCGGAACCA#","Polyribonucleotide (RNA)","","","","","","","","","","","","","","","","","","TYMV PSEUDOKNOT","","Turnip yellow mosaic virus","12154","","","","","","","","","","","","","","","","","","","","90.0","90.0","90.0","P 1","1.0","1.0","1.0","1","","","","","","","","","","","400.0","VARIAN","UNITY+","","","","","","","","","X-PLOR","","TORSION ANGLE DYNAMICS/ SIMULATED ANNEALING","REFINEMENT DETAILS CAN BE FOUND IN THE JRNL CITATION ABOVE.","LEAST RESTRAINT VIOLATIONS","140","24","","","","","",""', '']
        self.assertTrue(line, result)

    def test_load_first_few_pdbs_into_db(self):
        """Make sure there are no errors when loading into the db"""
        temp = self.loader.pdbs
        self.loader.pdbs = self.loader.pdbs[:3]
        status = self.loader.update_pdb_info()
        self.loader.pdbs = temp
        self.assertTrue(status)


if __name__ == '__main__':
    unittest.main()