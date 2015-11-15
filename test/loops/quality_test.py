from test import StageTest
from nose import SkipTest

from pymotifs.loops.quality import Loader


class QueryingTest(StageTest):
    loader_class = Loader

    def test_can_detect_if_has_data(self):
        self.assertTrue(self.loader.has_data('1GID'))

    def test_can_detect_if_has_no_data(self):
        self.assertFalse(self.loader.has_data('0GID'))

    def test_can_get_current_release(self):
        self.assertEquals('0.1', self.loader.current_id())


class LoopsTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(LoopsTest, self).setUp()
        self.loops = self.loader.loops('4V9Q')

    def test_can_load_all_loops(self):
        raise SkipTest()

    def test_builds_the_correct_loop(self):
        raise SkipTest()


class ComplementaryTest(StageTest):
    loader_class = Loader

    def test_knows_HL_not_complementary(self):
        loop = {'type': 'HL'}
        self.assertFalse(self.loader.is_complementary(loop))

    def test_knows_J3_not_complementary(self):
        loop = {'type': 'J3'}
        self.assertFalse(self.loader.is_complementary(loop))

    def test_knows_non_complemntary_IL(self):
        loop = {'type': 'IL', 'sequence': 'ACCAG,CGUGU'}
        self.assertFalse(self.loader.is_complementary(loop))

    def test_knows_an_even_length_is_complementary(self):
        loop = {'type': 'IL', 'sequence': 'ACGGUU,AACCGU'}
        self.assertTrue(self.loader.is_complementary(loop))

    def test_knows_an_odd_length_is_complementary(self):
        loop = {'type': 'IL', 'sequence': 'ACGUU,AACGU'}
        self.assertTrue(self.loader.is_complementary(loop))

    def test_knows_1_letter_not_complementary(self):
        loop = {'type': 'IL', 'sequence': 'A'}
        self.assertFalse(self.loader.is_complementary(loop))

    def test_knows_2_letter_not_complementary(self):
        loop = {'type': 'IL', 'sequence': 'AC'}
        self.assertFalse(self.loader.is_complementary(loop))

    def test_knows_2_letter_is_complementary(self):
        loop = {'type': 'IL', 'sequence': 'AU'}
        self.assertTrue(self.loader.is_complementary(loop))

    def test_knows_observed_are_complementary(self):
        known = ['CGU,AUG', 'CUA,UAG', 'GAGC,GCUC', 'GGC,GCC', 'AGC,GCU',
                 'GUG,CGC', 'UCA,UGG', 'CAG,CUG', 'AGA,UCU', 'UCC,GGA',
                 'GAC,GUC', 'GUG,UGC', 'GAU,AUC', 'GGA,UCC', 'GUG,CAC',
                 'CUG,CAG', 'CUGUC,GACAG', 'GGU,ACC', 'GUC,GAC', 'CCG,CGG',
                 'GAA,UUC', 'ACC,GGU', 'CUC,GGG', 'GUC,GGC', 'CGCG,CGUG',
                 'GGG,CUC', 'CGC,GUG', 'UGG,CUG', 'GCC,GGC', 'AUG,CGU',
                 'CUGG,CUGG', 'CGGU,AUCG', 'CGG,CUG', 'CUCG,CGAG', 'CUG,CGG',
                 'GUAC,GUGC', 'GUA,UGC', 'CAC,GUG', 'GGG,CCC', 'UGC,GCA',
                 'CCC,GGG', 'UAG,CUA', 'CUC,GAG', 'GCU,GGC', 'GGU,AUC',
                 'UUA,UAA', 'CGC,GCG', 'CGG,CCG', 'GCG,CGU', 'CGGA,UCCG',
                 'GAG,CUC', 'GGC,GCU', 'UGG,CCA', 'GGG,CCU', 'CAU,GUG',
                 'CCA,UGG', 'CGGG,CCUG', 'GCG,CGC', 'UCU,AGG', 'UGU,GCG',
                 'UAC,GUG', 'GGG,UCC', 'UGU,AUG', 'UCG,CGA', 'UCGC,GUGA',
                 'CGA,UUG', 'AUC,GGU', 'UUC,GGA', 'GACG,CGUC', 'GUG,CAU',
                 'GUAG,UUGC', 'UAG,UUG', 'AGA,UUU', 'GUAGA,UUUGC', 'UUG,CAA',
                 'GUG,CGU', 'GGU,GUC', 'GGC,GUC', 'GGG,CUU', 'CUG,UGG',
                 'AUCG,CGGU', 'GGA,UUU', 'GGCCGCG,CGUGGUC', 'CUCC,GGAG',
                 'CCCG,CGGG', 'CGCC,GGCG', 'CGGGC,GCCCG', 'GGCCGC,GUGGUC',
                 'GGGU,GCCU', 'CCCA,UGGG', 'CAA,UUG', 'UUCA,UGGA',
                 'AGU,GCU', 'UGG,UUA', 'AGC,GUU', 'AUG,UGU', 'GUAU,GUAC',
                 'CGGG,CUUG', 'GUUU,AAAC', 'UUG,UGA']
        fn = self.loader.is_complementary
        for seq in known:
            self.assertEquals(seq, fn({'type': 'IL', 'sequence': seq}))


class ModificationsTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(ModificationsTest, self).setUp()
        self.loops = {}
        for loop in self.loader.loops('1C2W'):
            self.loops[loop['id']] = loop

    def test_can_detect_a_modified_loop(self):
        raise SkipTest()
        val = self.loader.modifications(self.loops['HL_1C2W_038'])
        self.assertEquals(['PHI'], val)

    def test_can_detect_multiple_modifications(self):
        raise SkipTest()
        val = self.loader.modifications(self.loops['HL_1C2W_045'])
        self.assertEquals(['PHI'], val)

    def test_knows_unmodified_loop(self):
        raise SkipTest()
        val = self.loader.modifications(self.loops['HL_1C2W_044'])
        self.assertEquals([], val)


class IncompleteNucleotidesTest(StageTest):
    pass
