from test import StageTest

from pymotifs.utils import releases as rel


class ComputingNextReleaseIdTest(StageTest):
    def test_can_compute_the_next_id(self):
        self.assertEquals('1.73', rel.next_id('1.72'))

    def test_computes_next_major_id(self):
        self.assertEquals('2.0', rel.next_id('1.72', mode='major'))

    def test_computes_next_minor_id(self):
        self.assertEquals('1.73', rel.next_id('1.72', mode='minor'))

    def test_computes_next_no_change(self):
        self.assertEquals('1.72', rel.next_id('1.72', mode='none'))

    def test_complains_for_unknown_mode(self):
        self.assertRaises(rel.UnknownReleaseMode, rel.next_id, '1.2', mode='b')

    def test_complains_for_badly_formatted_release_id(self):
        self.assertRaises(rel.BadlyFormattedRelease, rel.next_id, '1')
        self.assertRaises(rel.BadlyFormattedRelease, rel.next_id, '1.')
        self.assertRaises(rel.BadlyFormattedRelease, rel.next_id, 'a.')
        self.assertRaises(rel.BadlyFormattedRelease, rel.next_id, '1.a')
