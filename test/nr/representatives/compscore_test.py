import pytest

from test import StageTest

from pymotifs.nr.representatives import CompScore


class DataTest(StageTest):
    loader_class = CompScore

    def data(self, *ids):
        return self.loader.load_quality([{'id': i for i in ids}])

    def quality(self, *args):
        return self.data(*args)['quality']

    def test_computes_correct_average_rsr(self):
        assert self.quality('1S72|1|0')['average_rsr'] == 0.13
        # assert self.quality('4v9f|0|1')['average_rsr'] == 0.14
        # assert self.quality('4v9f|9|2')['average_rsr'] == 0.13
        # assert self.quality('4V7M|32|DB')['average_rsr'] == 0.23

    def test_computes_correct_percent_clash(self):
        assert self.quality('1S72|1|0')['percent_clash'] == 0.06
        # assert self.quality('4v9f|0|1')['percent_clash'] == 0.07
        # assert self.quality('4v9f|9|2')['percent_clash'] == 0.15
        # assert self.quality('4V7M|32|DB')['percent_clash'] == 3.84

        assert 'percent_clash' in self.quality('1S72|1|0')['has']
        # assert 'percent_clash' in self.quality('4V7M|32|DB')['has']
        # assert 'percent_clash' in self.quality('4v9f|0|1')['has']
        # assert 'percent_clash' in self.quality('4v9f|9|2')['has']

    def test_computes_correct_average_rscc(self):
        assert self.quality('1S72|1|0')['average_rscc'] == -1 * (0.035 - 1)
        # assert self.quality('4v9f|0|1')['average_rscc'] == -1 * (0.033 - 1)
        # assert self.quality('4v9f|9|2')['average_rscc'] == -1 * (0.043 - 1)
        # assert self.quality('4V7M|32|DB')['average_rscc'] == -1 * (0.277 - 1)

    def test_uses_correct_rfree(self):
        assert self.quality('1S72|1|0')['average_rscc'] == 0.22
        # assert self.quality('4V7M|32|DB')['average_rscc'] == 0.27
        # assert self.quality('4v9f|0|1')['average_rscc'] == 0.21
        # assert self.quality('4v9f|9|2')['average_rscc'] == 0.21

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
