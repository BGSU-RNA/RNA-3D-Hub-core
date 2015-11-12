from test import StageTest

from pymotifs.utils.correspondence import Helper


class AlignedChainsTest(StageTest):
    loader_class = Helper

    def setUp(self):
        super(AlignedChainsTest, self).setUp()
        self.alignments = self.loader.aligned_chains(['1GID', '3T4B'])

    def test_can_load_all_alignments_from_the_pdb(self):
        self.assertEquals([60, 61, 172], sorted(self.alignments.keys()))

    def test_can_load_all_alignments_to_each(self):
        self.assertEquals({172: False, 61: True}, self.alignments[60])
