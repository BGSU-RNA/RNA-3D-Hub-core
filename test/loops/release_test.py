from test import StageTest

from pymotifs.utils.releases import Release
from pymotifs.utils.releases import UnknownReleaseType
from pymotifs.utils.releases import UnknownReleaseMode


class GettingReleaseIdTest(StageTest):
    loader_class = Release

    def test_can_get_current_loop_release(self):
        self.assertEquals('1.72', self.loader.current('loop'))

    def test_can_get_current_motif_release(self):
        self.assertEquals('1.16', self.loader.current('motif'))

    def test_can_get_current_nr_release(self):
        self.assertEquals('1.70', self.loader.current('nr'))

    def test_complains_for_unknown_release_type(self):
        self.assertRaises(UnknownReleaseType, self.loader.current, 'bob')


class ComputingNextReleaseIdTest(StageTest):
    loader_class = Release

    def test_can_compute_the_next_id(self):
        self.assertEquals('1.73', self.loader.next('1.72'))

    def test_computes_next_major_id(self):
        self.assertEquals('2.0', self.loader.next('1.72', mode='major'))

    def test_computes_next_minor_id(self):
        self.assertEquals('1.73', self.loader.next('1.72', mode='minor'))

    def test_computes_next_no_change(self):
        self.assertEquals('1.72', self.loader.next('1.72', mode='none'))

    def test_complains_for_unknown_mode(self):
        self.assertRaises(UnknownReleaseMode, self.loader.next, '1.72',
                          mode='bob')
