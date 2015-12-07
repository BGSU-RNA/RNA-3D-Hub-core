from test import StageTest
from nose import SkipTest

from pymotifs.nr.parent_counts import Loader


class KnownNamesTest(StageTest):
    loader_class = Loader

    def test_gets_all_given_resolution(self):
        raise SkipTest()
