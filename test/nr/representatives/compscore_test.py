import pytest

from numpy.testing import assert_almost_equal

from test import StageTest

from pymotifs.nr.representatives import CompScore


class DataTest(StageTest):
    loader_class = CompScore

    def data(self, *ids):
        return self.loader.load_quality([{'id': i for i in ids}])

    def quality(self, *args):
        return self.data(*args)[0]['quality']

    def test_computes_correct_average_rsr(self):
        assert_almost_equal(self.quality('1S72|1|0')['average_rsr'], 0.15, decimal=2)
        assert_almost_equal(self.quality('4V9F|0|1')['average_rsr'], 0.14, decimal=2)
        assert_almost_equal(self.quality('4V9F|9|2')['average_rsr'], 0.13, decimal=2)
        assert_almost_equal(self.quality('4V7M|32|DB')['average_rsr'], 0.23, decimal=2)

    def test_computes_correct_percent_clash(self):
        assert_almost_equal(self.quality('1S72|1|0')['percent_clash'], 0.056, decimal=2)
        assert_almost_equal(self.quality('4v9f|0|1')['percent_clash'], 0.07, decimal=2)
        assert_almost_equal(self.quality('4v9f|9|2')['percent_clash'], 0.15, decimal=2)
        assert_almost_equal(self.quality('4V7M|32|DB')['percent_clash'], 3.84, decimal=2)

        assert 'percent_clash' in self.quality('1S72|1|0')['has']
        assert 'percent_clash' in self.quality('4V7M|32|DB')['has']
        assert 'percent_clash' in self.quality('4v9f|0|1')['has']
        assert 'percent_clash' in self.quality('4v9f|9|2')['has']

    def test_computes_correct_average_rscc(self):
        assert_almost_equal(self.quality('1S72|1|0')['average_rscc'], 0.961, decimal=2)
        assert self.quality('4v9f|0|1')['average_rscc'] == -1 * (0.033 - 1)
        assert self.quality('4v9f|9|2')['average_rscc'] == -1 * (0.043 - 1)
        assert self.quality('4V7M|32|DB')['average_rscc'] == -1 * (0.277 - 1)

    def test_uses_correct_rfree(self):
        assert self.quality('1S72|1|0')['rfree'] == 0.22
        assert self.quality('4V7M|32|DB')['rfree'] == 0.27
        assert self.quality('4v9f|0|1')['rfree'] == 0.21
        assert self.quality('4v9f|9|2')['rfree'] == 0.21


    @pytest.mark.skip()
    def test_complains_about_missing_rfree(self):
        with pytest.raises(Exception):
            self.data('157D|0|A')

    def test_computes_correct_compscore(self):
        assert self.loader.compscore(self.data('1S72|1|0')) == 127
        # assert self.loader.compscore(self.data('4v9f|0|1')) == 125
        # assert self.loader.compscore(self.data('4v9f|9|2')) == 126
        # assert self.loader.compscore(self.data('4V7M|32|DB')) == 127


class SelectingRepresentativeTest(StageTest):
    loader_class = CompScore

    @pytest.mark.skip()
    def test_selects_correct_representative(self):
        members = [
            self.data('4v9f|9|2'),
            self.data('1S72|0|1'),
            self.data('4v9f|0|1'),
            self.data('4V7M|32|DB'),
        ]
        assert self.representative(members) == self.data('4v9f', 0, '1')


class SpecificExamples(StageTest):
    loader_class = CompScore

    def compscore(self, ife_id):
        info = {'id': ife_id}
        members = self.loader.load_quality([info])
        return self.loader.compscore(members[0])

    def test_157D_A(self):
        assert self.compscore('157D|1|A+157D|1|B') == 216.0

    def test_1CGM(self):
        assert self.compscore('1CGM|1|I') == 2223.33
