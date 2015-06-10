from test import StageTest

from pymotifs.utils.correspondence import Helper


class AlignedChainsTest(StageTest):
    loader_class = Helper

    def setUp(self):
        super(AlignedChainsTest, self).setUp()
        self.alignments = self.loader.aligned_chains(['4KJ2', '4KJ4'])

    def test_can_load_all_alignments_from_the_pdb(self):
        self.assertEquals([24284, 24230], self.alignments.keys())

    def test_can_load_all_alignments_to_each(self):
        self.assertEquals({24284: 1}, self.alignments[24230])
