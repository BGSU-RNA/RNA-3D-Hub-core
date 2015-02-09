import unittest as ut

from pymotifs.utils.pdb import CustomReportHelper


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
