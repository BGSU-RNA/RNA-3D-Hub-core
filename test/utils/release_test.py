from test import StageTest

from pymotifs.utils.releases import Release
from pymotifs.utils.releases import UnknownReleaseType
from pymotifs.utils.releases import UnknownReleaseMode
from pymotifs.utils.releases import BadlyFormattedRelease


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


class GettingWithNoCurrentReleaseTest(StageTest):
    loader_class = Release

    def setUp(self):
        super(GettingWithNoCurrentReleaseTest, self).setUp()
        with self.loader.session() as session:
            session.execute('CREATE TABLE ml_releases_tmp AS (SELECT * FROM ml_releases);')
            session.execute('delete from ml_releases')

    def tearDown(self):
        with self.loader.session() as session:
            session.execute('insert into ml_releases select * from ml_releases_tmp ;')
            session.execute('drop table ml_releases_tmp')

    def test_sets_missing_release_to_zero(self):
        self.assertEquals('0.0', self.loader.current('motif'))


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

    def test_complains_for_badly_formatted_release_id(self):
        self.assertRaises(BadlyFormattedRelease, self.loader.next, '1')
        self.assertRaises(BadlyFormattedRelease, self.loader.next, '1.')
        self.assertRaises(BadlyFormattedRelease, self.loader.next, 'a.')
        self.assertRaises(BadlyFormattedRelease, self.loader.next, '1.a')
