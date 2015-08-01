from nose import SkipTest

from pymotifs.core import Stage

from test import StageTest as Base
from test import config


class SomeStage(Stage):
    def is_missing(self, data, **kwargs):
        return data


class StageNameTest(Base):
    def test_it_can_generate_a_name(self):
        val = SomeStage(config, None).name
        ans = 'test.core.stage_test'
        self.assertEqual(ans, val)


class LoggerTest(Base):
    def test_it_creates_a_logger(self):
        val = SomeStage(config, None)
        self.assertTrue(val.logger)

    def test_it_assigns_the_loggers_name(self):
        val = SomeStage(config, None)
        self.assertEqual('test.core.stage_test', val.logger.name)


class MustRecomputeTest(Base):
    def test_must_recompute_detects_if_given_recompute(self):
        conf = dict(config)
        conf.update({'test.core.stage_test': {'recompute': False}})
        stage = SomeStage(conf, None)
        val = stage.must_recompute(None, recalculate=True)
        self.assertTrue(val)

    def test_must_recompute_detects_config_has_recompute(self):
        conf = dict(config)
        conf.update({'test.core.stage_test': {'recompute': True}})
        stage = SomeStage(conf, None)
        val = stage.must_recompute(None, recalculate=False)
        self.assertTrue(val)

    def test_will_not_recompute_otherwise(self):
        conf = dict(config)
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
        stage = SomeStage(config, None)
        val = stage.should_process(None, recalculate=True)
        self.assertTrue(val)

    def test_it_will_reprocess_if_missing_data(self):
        stage = SomeStage(config, None)
        val = stage.should_process(True)
        self.assertTrue(val)

    def test_will_reprocess_if_config_has_recompute(self):
        conf = dict(config)
        conf.update({'test.core.stage_test': {'recompute': True}})
        stage = SomeStage(conf, None)
        val = stage.should_process(None)
        self.assertTrue(val)

    def test_it_will_not_reprocess_otherwise(self):
        stage = SomeStage(config, None)
        val = stage.should_process(None)
        self.assertFalse(val)


class ProcessingTests(Base):
    def test_it_will_convert_intput_to_upper(self):
        stage = SomeStage(config, None)
        val = stage.to_process(['abc'])
        ans = ['ABC']
        self.assertEqual(val, ans)
