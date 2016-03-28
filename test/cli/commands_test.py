import random

import pytest

from unittest import TestCase

from pymotifs.cli import commands
from pymotifs.cli.commands import Runnable

from test import CONFIG

commands.SHOULD_RUN = False


class RunnableTest(TestCase):

    def do(self, config=None, options=None):
        engine = None
        config = config or CONFIG
        options = options or {}
        runner = Runnable('units.info', ['1GID'], config, engine, options)
        return runner, runner()

    def names(self, runner, name):
        stages = runner.dispatcher.stages(name, build=True)
        return [s.name.replace('pymotifs.', '') for s in stages]

    def recalc(self, name):
        pass

    @pytest.mark.skip()
    def test_will_run_all_stages(self):
        assert False

    def test_will_set_the_seed_if_given(self):
        # Kinda tricky to test if the seed is set directly, instead we test if
        # the in generated is the same as what happens if you set the seed
        # manually. This will break if the generator method changes in the
        # future.
        self.do(options={'seed': 1})
        assert random.randint(0, 10) == 1

    def test_it_will_not_skip_stages_unrequested(self):
        runner, _ = self.do(options={'skip_stage': []})
        assert runner.skipped_stages == []
        assert self.names(runner, 'units.info') == \
            ['download', 'pdbs.info', 'units.info']
        runner, _ = self.do()
        assert runner.skipped_stages == []

    def test_will_skip_marked_stages(self):
        runner, _ = self.do(options={'skip_stage': ['download']})
        assert runner.skipped_stages == ['download']
        assert self.names(runner, 'units.info') == ['pdbs.info', 'units.info']

    def test_it_will_skip_current_stage_if_asked(self):
        runner, _ = self.do(options={'skip_stage': ['.']})
        assert runner.skipped_stages == ['units.info']
        assert 'units.info' not in self.names(runner, 'units.info')

    def test_will_recalculate_all_stages_given_star(self):
        runner, _ = self.do(options={'recalculate': ['*']})
        assert runner.recalculate_stages is True
        runner, _ = self.do(options={'recalculate': ['.', '*']})
        assert runner.recalculate_stages is True

    def test_will_recalculate_current_stage_given_dot(self):
        runner, _ = self.do(options={'recalculate': ['.']})
        assert runner.recalculate_stages == ['units.info']

    def test_it_will_not_recalculate_any_stages_if_not_set(self):
        runner, _ = self.do()
        assert runner.recalculate_stages == []

    def test_will_recalculate_given_stages(self):
        runner, _ = self.do(options={'recalculate': ['.', 'download']})
        assert runner.recalculate_stages == ['units.info', 'download']

    @pytest.mark.skip()
    def test_will_send_failing_email_if_stage_fails(self):
        pass

    @pytest.mark.skip()
    def test_will_send_suceeded_email_if_stage_works(self):
        pass

    @pytest.mark.skip()
    def test_will_not_send_mail_if_not_set(self):
        pass
