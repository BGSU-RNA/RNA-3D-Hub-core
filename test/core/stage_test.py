import os

import pytest

from pymotifs.core.stages import Stage
from pymotifs.core import Skip

from test import StageTest as Base
from test import CONFIG


class SomeStage(Stage):
    def is_missing(self, data, **kwargs):
        return data

    def store(self, data, **kwargs):
        pass

    def mark_processed(self, data, **kwargs):
        pass


class RecomputingTest(Base):
    def test_defaults_to_not_recomputing(self):
        stage = SomeStage(CONFIG, None)
        val = stage.must_recompute(None)
        self.assertFalse(val)

    def test_must_recompute_detects_if_given_recompute(self):
        conf = dict(CONFIG)
        conf.update({'test.core.stage_test': {'recompute': False}})
        stage = SomeStage(conf, None)
        val = stage.must_recompute(None, recalculate=True)
        self.assertTrue(val)

    def test_will_not_recompute_if_given_false(self):
        stage = SomeStage(CONFIG, None)
        val = stage.must_recompute(None, recalculate=False)
        self.assertFalse(val)

    def test_must_recompute_detects_config_has_recompute(self):
        conf = dict(CONFIG)
        conf.update({'test.core.stage_test': {'recompute': True}})
        stage = SomeStage(conf, None)
        val = stage.must_recompute(None, recalculate=False)
        self.assertTrue(val)

    def test_will_not_recompute_otherwise(self):
        conf = dict(CONFIG)
        conf.update({'test.core.stage_test': {'recompute': False}})
        stage = SomeStage(conf, None)
        val = stage.must_recompute(None, recalculate=False)
        self.assertFalse(val)

    def test_will_not_recompute_if_in_skip(self):
        stage = SomeStage(CONFIG, None, skip_pdbs=['0GGG'])
        self.assertRaises(Skip, stage.should_process, '0GGG')

    def test_will_recompute_if_given_true(self):
        stage = SomeStage(CONFIG, None)
        val = stage.must_recompute(None, recalculate=True)
        self.assertTrue(val)

    def test_will_recompute_if_given_list_with_name(self):
        stage = SomeStage(CONFIG, None)
        recalc = ['test.core.stage_test', 'units.info']
        val = stage.must_recompute(None, recalculate=recalc)
        self.assertTrue(val)

    def test_will_not_recompute_if_name_not_in_list(self):
        stage = SomeStage(CONFIG, None)
        recalc = ['test.core.bob', 'units.info']
        val = stage.must_recompute(None, recalculate=recalc)
        self.assertFalse(val)

    def test_will_recompute_if_given_set_with_name(self):
        stage = SomeStage(CONFIG, None)
        recalc = set(['test.core.stage_test', 'units.info'])
        val = stage.must_recompute(None, recalculate=recalc)
        self.assertTrue(val)

    def test_will_not_recompute_if_name_not_in_set(self):
        stage = SomeStage(CONFIG, None)
        recalc = set(['test.core.bob', 'units.info'])
        val = stage.must_recompute(None, recalculate=recalc)
        self.assertFalse(val)

    def test_will_recompute_if_given_tuple_with_name(self):
        stage = SomeStage(CONFIG, None)
        recalc = ('test.core.stage_test', 'units.info')
        val = stage.must_recompute(None, recalculate=recalc)
        self.assertTrue(val)

    def test_will_not_recompute_if_name_not_in_tuple(self):
        stage = SomeStage(CONFIG, None)
        recalc = ('test.core.bob', 'units.info')
        val = stage.must_recompute(None, recalculate=recalc)
        self.assertFalse(val)


class BeenLongEnoughTest(Base):
    @pytest.mark.skip()
    def test_knows_if_it_is_too_long(self):
        pass

    @pytest.mark.skip()
    def test_knows_if_not_too_long(self):
        pass


class ShouldProcessTest(Base):
    @pytest.mark.skip()
    def test_will_reprocess_if_too_long(self):
        pass

    def test_it_will_reprocess_if_given_recalculate(self):
        stage = SomeStage(CONFIG, None)
        val = stage.should_process(None, recalculate=True)
        self.assertTrue(val)

    def test_it_will_reprocess_if_missing_data(self):
        stage = SomeStage(CONFIG, None)
        val = stage.should_process(True)
        self.assertTrue(val)

    def test_will_reprocess_if_config_has_recompute(self):
        conf = dict(CONFIG)
        conf.update({'test.core.stage_test': {'recompute': True}})
        stage = SomeStage(conf, None)
        val = stage.should_process(None)
        self.assertTrue(val)

    def test_it_will_not_reprocess_otherwise(self):
        stage = SomeStage(CONFIG, None)
        val = stage.should_process(None)
        self.assertFalse(val)


class ProcessingTests(Base):
    def test_it_will_convert_intput_to_upper(self):
        stage = SomeStage(CONFIG, None)
        val = stage.to_process(['abc'])
        ans = ['ABC']
        self.assertEqual(val, ans)

    def test_it_will_convert_to_str(self):
        stage = SomeStage(CONFIG, None)
        val = stage.to_process([u'abc'])
        self.assertEqual(['ABC'], val)
        self.assertTrue(isinstance(val[0], str))

    def test_will_return_processed_input(self):
        stage = SomeStage(CONFIG, None)
        val = stage(['A', '', 'B'])
        self.assertEqual(val, ['A', 'B'])


class CachingTest(Base):
    loader_class = SomeStage

    def tearDown(self):
        filename = self.loader.cache_filename('example')
        if os.path.exists(filename):
            os.remove(filename)

    def test_it_can_cache_data(self):
        ans = [1]
        self.loader.cache('example', ans)
        self.assertEquals(ans, self.loader.cached('example'))

    def test_it_gives_none_for_no_cached_data(self):
        self.assertEquals(None, self.loader.cached('bob'))

    def test_it_can_remove_cached_data(self):
        ans = [1]
        self.loader.cache('example', ans)
        self.assertEquals(ans, self.loader.cached('example', remove=True))
        self.assertEquals(None, self.loader.cached('example', remove=True))
