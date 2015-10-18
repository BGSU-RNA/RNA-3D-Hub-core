from test import StageTest

from datetime import date
from datetime import timedelta

from pymotifs.models import NrReleases
from pymotifs.utils.releases import Release
from pymotifs.utils.releases import UnknownReleaseType
from pymotifs.utils.releases import UnknownReleaseMode
from pymotifs.utils.releases import BadlyFormattedRelease
from pymotifs.utils.releases import ImpossibleRelease


class GettingReleaseIdTest(StageTest):
    loader_class = Release

    def tearDown(self):
        with self.loader.session() as session:
            session.query(NrReleases).delete()

    def test_can_get_current_loop_release(self):
        self.assertEquals('0.0', self.loader.current('loop'))

    def test_can_get_current_motif_release(self):
        self.assertEquals('0.0', self.loader.current('motif'))

    def test_can_get_current_nr_release(self):
        self.assertEquals('0.0', self.loader.current('nr'))

    def test_can_get_current_nr_release_with_existing_one(self):
        today = date.today()
        yesterday = today - timedelta(1)
        with self.loader.session() as session:
            session.add(NrReleases(nr_release_id='0.1', date=yesterday))
            session.add(NrReleases(nr_release_id='1.1', date=today))

        self.assertEquals('1.1', self.loader.current('nr'))

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
            session.execute('insert into ml_releases select * from ml_releases_tmp;')
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


class GettingPreviousReleaseId(StageTest):
    loader_class = Release

    def test_can_get_previous_release(self):
        self.assertEquals('1.72', self.loader.previous('1.73'))

    def test_can_compute_previous_major_release(self):
        self.assertEquals('1.0', self.loader.previous('2.54', mode='major'))

    def test_computes_previous_minor_id(self):
        self.assertEquals('3.2', self.loader.previous('3.3', mode='minor'))

    def test_does_not_wrap_major_below_zero(self):
        self.assertEquals('0.1', self.loader.previous('0.54', mode='major'))

    def test_does_not_wrap_minor_below_zero(self):
        self.assertEquals('1.0', self.loader.previous('1.0', mode='minor'))

    def test_will_give_none_if_requested_for_minor(self):
        self.assertEquals(None, self.loader.previous('1.0', mode='minor',
                                                     bad_id='none'))

    def test_will_give_none_if_requested_for_major(self):
        self.assertEquals(None, self.loader.previous('0.0', mode='major',
                                                     bad_id='none'))

    def test_will_raise_if_requested_for_bad_minor(self):
        self.assertRaises(ImpossibleRelease, self.loader.previous,
                          '1.0', mode='minor', bad_id='raise')

    def test_will_raise_if_requested_for_bad_major(self):
        self.assertRaises(ImpossibleRelease, self.loader.previous,
                          '0.0', mode='major', bad_id='raise')

    # def test_can_get_previous_from_database(self):
    #     self.assertEquals('0.110', self.loader.previous('1.0', mode='lookup',
    #                                                     name='nr'))

    def test_it_gives_none_if_asked_for_previous_of_unknown(self):
        self.assertEquals(None, self.loader.previous('3.0', mode='lookup',
                                                     name='nr'))
