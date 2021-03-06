from test import StageTest
import hashlib

from pymotifs.exp_seq.info import Loader


class InfoTest(StageTest):
    loader_class = Loader

    def test_knows_if_it_has_data(self):
        self.assertTrue(self.loader.has_data('UUUUCU'))

    def test_knows_if_it_does_not_have_data(self):
        self.assertFalse(self.loader.has_data('UUUUGU'))


class ToProcessTest(StageTest):
    loader_class = Loader

    def test_can_get_all_unique_sequences(self):
        self.assertEquals(2, len(self.loader.to_process(['1S72'])))

    def test_loads_only_unique_sequences(self):
        val = self.loader.to_process(['157D', '124D'])
        ans = ['CAUGUGAC', 'CGCGAAUUAGCG']
        self.assertEquals(ans, val)


class NormalizationTest(StageTest):
    loader_class = Loader

    def test_can_normalize_a_sequence(self):
        self.assertEquals('ACGUNN', self.loader.normalize('ACIUNN'))

    def test_does_not_modify_a_good_sequence(self):
        self.assertEquals('ACGUNN', self.loader.normalize('ACGUNN'))

    def test_turns_X_into_N(self):
        self.assertEquals('ACGUNN', self.loader.normalize('ACGUNX'))

    def test_turns_I_into_G(self):
        self.assertEquals('ACGUNN', self.loader.normalize('ACIUNN'))

    def test_returns_none_if_cannot_normalize_all_characters(self):
        self.assertEquals(None, self.loader.normalize('ACGBN'))

    def test_will_remove_last_char_if_cannot_normalize(self):
        self.assertEquals('ACGUNN', self.loader.normalize('ACGUNXT'))


class DataTest(StageTest):
    loader_class = Loader

    def test_builds_correctly_normalized_data(self):
        seq = 'ACIUNXT'
        norm = 'ACGUNN'
        val = self.loader.data(seq)
        ans = {
            'sequence': seq,
            'md5': hashlib.md5(seq).hexdigest(),
            'normalized': norm,
            'length': len(seq),
            'normalized_length': len(norm),
            'was_normalized': True
        }
        self.assertEquals(ans, val)

    def test_indicates_not_normalized(self):
        seq = 'ACI?NXT'
        val = self.loader.data(seq)
        ans = {
            'sequence': seq,
            'md5': hashlib.md5(seq).hexdigest(),
            'normalized': None,
            'length': len(seq),
            'normalized_length': 0,
            'was_normalized': False
        }
        self.assertEquals(ans, val)
