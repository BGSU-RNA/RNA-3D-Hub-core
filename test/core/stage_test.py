from nose import SkipTest

from pymotifs.core import Stage

from test import StageTest as Base


class SomeStage(Stage):
    name = 'bob'

    def is_missing(self, data, **kwargs):
        return data


class NamelessStage(Stage):
    pass


class StageNameTest(Base):
    def test_it_uses_the_given_name(self):
        val = SomeStage({}, None).name
        ans = 'bob'
        self.assertEqual(ans, val)

    def test_it_can_generate_a_name(self):
        val = NamelessStage({}, None).name
        ans = 'test.core.stage_test'
        self.assertEqual(ans, val)


def LoggerTest(Base):
    def test_it_creates_a_logger(self):
        val = SomeStage({}, None)
        self.assertTrue(val.logger)

    def test_it_assigns_the_loggers_name(self):
        val = SomeStage({}, None)
        self.assertTrue('bob', val.logger.name)


class MustRecomputeTest(Base):
    def test_must_recompute_detects_if_given_recompute(self):
        stage = SomeStage({'bob': {'recompute': False}}, None)
        val = stage.must_recompute(None, recalculate=True)
        self.assertTrue(val)

    def test_must_recompute_detects_config_has_recompute(self):
        stage = SomeStage({'bob': {'recompute': True}}, None)
        val = stage.must_recompute(None, recalculate=False)
        self.assertTrue(val)

    def test_will_not_recompute_otherwise(self):
        stage = SomeStage({'bob': {'recompute': False}}, None)
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
        stage = SomeStage({}, None)
        val = stage.should_process(None, recalculate=True)
        self.assertTrue(val)

    def test_it_will_reprocess_if_missing_data(self):
        stage = SomeStage({}, None)
        val = stage.should_process(True)
        self.assertTrue(val)

    def test_will_reprocess_if_config_has_recompute(self):
        stage = SomeStage({'bob': {'recompute': True}}, None)
        val = stage.should_process(None)
        self.assertTrue(val)

    def test_it_will_not_reprocess_otherwise(self):
        stage = SomeStage({}, None)
        val = stage.should_process(None)
        self.assertFalse(val)


class ProcessingTests(Base):
    def test_it_will_convert_intput_to_upper(self):
        stage = SomeStage({}, None)
        val = stage.to_process(['abc'])
        ans = ['ABC']
        self.assertEqual(val, ans)
