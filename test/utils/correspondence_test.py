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
        ans = {
            172L: {60L: False, 61L: False},
            60L: {172L: False, 61L: True},
            61L: {172L: False, 60L: True}
        }
        self.assertEquals(ans, self.alignments)

    def test_can_load_only_good_alignments(self):
        val = self.loader.aligned_chains(['1GID', '3T4B'], good=True)
        ans = {60L: {61L: True}, 61L: {60L: True}}
        self.assertEquals(ans, val)

    def test_can_load_only_bad_alignments(self):
        val = self.loader.aligned_chains(['1GID', '3T4B'], good=False)
        ans = {
            172L: {60L: False, 61L: False},
            60L: {172L: False},
            61L: {172L: False}
        }
        self.assertEquals(ans, val)
