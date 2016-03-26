from unittest import TestCase

import pytest

from pymotifs.cli import setup
from pymotifs.cli import introspect


def expand(name, value):
    return setup.expand_stage_pattern(name, value)


class ExpandStagePatternTest(TestCase):
    def test_expands_dot(self):
        assert expand('units.info', '.') == ['units.info']
        assert expand('units.info', ['.', '.']) == ['units.info', 'units.info']

    def test_expands_star(self):
        assert expand('units.info', '*') is True
        assert expand('units.info', ['.', '*', '.']) is True

    def test_gives_empty_if_nothing(self):
        assert setup.expand_stage_pattern('units.info', None) == []
        assert setup.expand_stage_pattern('units.info', []) == []

    def test_fails_given_unknown_stage(self):
        with pytest.raises(introspect.UnknownStageError):
            expand('units.info', ['bob'])

    def test_it_passes_names_through(self):
        assert expand('units.info', ['units.loader']) == ['units.loader']
        assert expand('units.info', ['pdbs.obsolete', 'units.quality']) == \
            ['pdbs.obsolete', 'units.quality']
