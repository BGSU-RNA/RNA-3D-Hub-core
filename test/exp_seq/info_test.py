from nose import SkipTest

from test import StageTest

from pymotifs.exp_seq.info import Loader


class InfoTest(StageTest):
    loader_class = Loader

    def test_can_load_all_sequences(self):
        raise SkipTest()
