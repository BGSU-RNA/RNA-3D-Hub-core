from test import StageTest
from test import config

from pymotifs.core import Base


class Simple(Base):
    pass


class BasicInstancePropertiesTest(StageTest):
    def test_it_creates_a_logger(self):
        val = Simple(config, None)
        self.assertTrue(val.logger)

    def test_it_assigns_the_loggers_name(self):
        val = Simple(config, None)
        self.assertEqual('test.core.base_test', val.logger.name)

    def test_name_defaults_to_module(self):
        val = Simple(config, None)
        self.assertEquals(val.name, 'test.core.base_test')
