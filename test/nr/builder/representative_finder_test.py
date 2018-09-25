import pytest

from pymotifs import core

from test import StageTest

from pymotifs.nr.builder import RepresentativeFinder


class RepresentativeFinderTest(StageTest):
    loader_class = RepresentativeFinder

    def test_can_get_object_of_method(self):
        assert self.loader.method('naive').method == 'naive'

    def test_will_fail_with_no_group(self):
        with pytest.raises(core.InvalidState):
            self.loader({})

    def test_will_fail_with_bad_name(self):
        with pytest.raises(core.InvalidState):
            self.loader({'id': 'a'}, method='impossible')

    @pytest.mark.skip()
    def test_can_pick_representative_using_method(self):
        pass
