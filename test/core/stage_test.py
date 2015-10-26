from nose import SkipTest

from pymotifs.core.stages import Stage

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
    def test_must_recompute_detects_if_given_recompute(self):
        conf = dict(CONFIG)
        conf.update({'test.core.stage_test': {'recompute': False}})
        stage = SomeStage(conf, None)
        val = stage.must_recompute(None, recalculate=True)
        self.assertTrue(val)

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


class BeenLongEnoughTest(Base):
    def test_knows_if_it_is_too_long(self):
        raise SkipTest()

    def test_knows_if_not_too_long(self):
        raise SkipTest()


class ShouldProcessTest(Base):
    def test_will_reprocess_if_too_long(self):
        raise SkipTest()

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

    def test_will_return_processed_input(self):
        stage = SomeStage(CONFIG, None)
        val = stage(['A', '', 'B'])
        self.assertEqual(val, ['A', 'B'])
