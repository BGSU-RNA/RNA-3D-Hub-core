import unittest as ut

from pymotifs.utils.pdb import CustomReportHelper
from pymotifs.chains.info import Loader


class CustomReportHelperTest(ut.TestCase):

    def setUp(self):
        self.helper = CustomReportHelper(fields=['structureId', 'ndbId'])

    def test_can_get_a_report_for_one_pdb(self):
        report = self.helper('2AW7')
        self.assertEquals(report, [{'structureId': '2AW7', 'ndbId': 'RR0125'}])

    def test_can_get_a_report_for_several_pdbs(self):
        report = self.helper(['157D', '2AW7'])
        val = [r['ndbId'] for r in report]
        self.assertEquals(2, len(report))
        self.assertEquals(['ARL048', 'RR0125'], val)

    def test_can_get_a_correct_report_for_1S72(self):
        self.helper = CustomReportHelper(fields=Loader.names.keys())
        report = self.helper(['1S72'])
        self.assertEquals(31, len(report))
