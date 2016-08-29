import os
import collections as coll
from unittest import TestCase

from pymotifs import config


class ConfigTest(TestCase):
    def setUp(self):
        self.conf = config.load("conf/bootstrap.json.txt")

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

    def test_it_keeps_known_keys(self):
        assert set(self.conf['locations'].keys()) == set([
            "2ds_destination",
            "cache",
            "fr3d_root",
            "interactions_gz",
            "log_dir",
            "loops_gz",
            "loops_mat_files",
            "loops_search_dir",
            "mlab_app",
            "releases_dir",
            '2ds_destination',
            'base',
            'cache',
            'fr3d_root',
            'log_dir',
            'loops_mat_files',
            'loops_search_dir',
            'releases_dir',
        ])

    def test_does_not_override_nested_values(self):
        val = self.conf['locations']['loops_gz']
        assert val == "/home/pipeline/hub-core/MotifAtlas/loops.gz"


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

    def will_merge_both_directions(self):
        a = {'a': {'a': {'c': 1}}}
        b = {'a': {'a': {'d': 1}, 'b': 3}}
        assert config.merge(a, b) == {'a': {'a': {'c': 1}, 'b': 3}}
        assert config.merge(b, a) == {'a': {'a': {'c': 1}, 'b': 3}}

    def test_will_preserve_a_default_dict(self):
        a = {'a': {'a': 1}}
        b = {'a': {'b': coll.defaultdict(int)}}
        ans = {'a': {'a': 1, 'b': coll.defaultdict(int)}}
        assert config.merge(a, b) == ans
        assert config.merge(b, a) == ans

    def test_can_merge_values_into_a_default_dict(self):
        def subdict(**kwargs):
            subdict = coll.defaultdict(int)
            for key, value in kwargs.items():
                subdict[key] = value
            return subdict

        a = {'a': {'a': 1, 'b': subdict()}}
        b = {'a': {'b': {'c': 5}}}
        ans = {'a': {'a': 1, 'b': subdict(c=5)}}
        assert config.merge(a, b) == ans
        assert config.merge(b, a) == ans
