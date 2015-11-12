import os
from unittest import TestCase

from pymotifs import config


class ConfigTest(TestCase):
    def setUp(self):
        self.conf = config.load("conf/test.json")

    def test_can_load_data(self):
        self.assertTrue('db' in self.conf)

    def test_uses_default_dict(self):
        val = self.conf['locations']['log_dir']
        ans = os.path.join(os.getcwd(), "MotifAtlas", "logs")
        self.assertEquals(ans, val)

    def test_adds_base_directory(self):
        val = self.conf['locations']['base']
        ans = os.getcwd()
        self.assertEquals(ans, val)

    def test_does_not_override_nested_values(self):
        val = self.conf['locations']['loops_mat_files']
        ans = "/home/pipeline/hub-core/MotifAtlas/PrecomputedData"
        self.assertEquals(ans, val)


class MergingTest(TestCase):
    def test_can_merge_two_dicts(self):
        a = {'c': 3}
        b = {'d': [1, 2, 3]}
        val = config.merge(a, b)
        ans = {'c': 3, 'd': [1, 2, 3]}
        self.assertEquals(ans, val)

    def test_can_recursively_merge(self):
        a = {'a': {'a': 1}, 'c': 3}
        b = {'a': {'b': 2}, 'd': [1, 2, 3]}
        val = config.merge(a, b)
        ans = {'a': {'a': 1, 'b': 2}, 'c': 3, 'd': [1, 2, 3]}
        self.assertEquals(ans, val)

    def test_converts_value_to_str(self):
        a = {'a': {'a': u'1'}}
        b = {'a': {'b': u'2'}}
        val = config.merge(a, b)
        self.assertTrue(isinstance(val['a']['b'], str))
