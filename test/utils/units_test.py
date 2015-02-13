from test import Session

from unittest import TestCase

from pymotifs.utils import units


class TranslatorTest(TestCase):

    def setUp(self):
        self.trans = units.Translator(Session)

    def test_can_translate_a_known_id(self):
        self.assertEquals(['1KOG|1|P|U|82'],
                          self.trans.translate('1KOG_AU_1_P_82_U_'))

    def test_can_translate_a_comma_seperated_string(self):
        ids = ['1KOG_AU_1_P_82_U_', '1KOG_AU_1_P_83_U_', '1KOG_AU_1_P_84_U_',
               '1KOG_AU_1_P_85_C_', '1KOG_AU_1_P_86_G_']
        ans = ['1KOG|1|P|U|82', '1KOG|1|P|U|83', '1KOG|1|P|U|84',
               '1KOG|1|P|C|85', '1KOG|1|P|G|86']
        self.assertEquals(ans, self.trans.translate(','.join(ids)))

    def test_can_translate_a_list_of_ids(self):
        ids = ['1KOG_AU_1_P_82_U_', '1KOG_AU_1_P_83_U_', '1KOG_AU_1_P_84_U_',
               '1KOG_AU_1_P_85_C_', '1KOG_AU_1_P_86_G_']
        ans = ['1KOG|1|P|U|82', '1KOG|1|P|U|83', '1KOG|1|P|U|84',
               '1KOG|1|P|C|85', '1KOG|1|P|G|86']
        self.assertEquals(ans, self.trans.translate(ids))

    def test_complains_if_missing_the_id(self):
        self.assertRaises(units.TranslationFailed, self.trans.translate, 'c')

    def test_complains_if_missing_one_id(self):
        ids = ['1KOG_AU_1_P_82_U_', '1KOG_AU_1_P_83_U', '1KOG_AU_1_P_84_U_',
               '1KOG_AU_1_P_85_C_', '1KOG_AU_1_P_86_G_']
        self.assertRaises(units.TranslationFailed, self.trans.translate, ids)
