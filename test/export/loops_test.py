from test import StageTest

from pymotifs.export.loops import Exporter


class GettingLoopsTest(StageTest):
    loader_class = Exporter

    def test_will_get_all_loops(self):
        data = self.loader.loops('1GID')
        assert len(data) == 22

    def test_creates_correct_data(self):
        data = self.loader.loops('1GID')
        assert data[0] == {
            'id': 'HL_1GID_001',
            'motif_id': None,
            'pdb': '1GID',
            'nts': '1GID|1|A|A|151,1GID|1|A|A|152,1GID|1|A|A|153,1GID|1|A|C|154,1GID|1|A|G|149,1GID|1|A|G|150',
            # 'motif_id': "HL_67042.17"
        }
