import unittest as ut

from models import Session
from correspondence.interactions import Loader


class InteractionLoaderTest(ut.TestCase):
    def setUp(self):
        self.loader = Loader(Session, {})

    def test_can_get_mapping(self):
        mapping = self.loader.mapping(38610)
        val = mapping['1FJG_AU_1_A_5_U_']
        ans = '3V22_BA1_1_A_5_U_'
        self.assertEqual(ans, val)
