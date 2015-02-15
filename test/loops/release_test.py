from test import StageTest
from nose import SkipTest


class GettingReleaseIdTest(StageTest):
    loader_class = StageTest

    def test_can_compute_the_next_id(self):
        raise SkipTest()

    def test_computes_next_major_id(self):
        raise SkipTest()

    def test_computes_next_minor_id(self):
        raise SkipTest()
