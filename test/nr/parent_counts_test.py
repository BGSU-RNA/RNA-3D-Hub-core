import pytest

from test import StageTest

from pymotifs.nr.parent_counts import Loader


class KnownNamesTest(StageTest):
    loader_class = Loader

    @pytest.mark.skip()
    def test_gets_all_given_resolution(self):
        pass
