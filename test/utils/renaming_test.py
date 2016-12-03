import unittest as ut

import pytest

import pymotifs.utils.renaming as rn


@pytest.mark.parametrize('test_input', [
    'a',
    '       a',
    '  a  ',
])
def test_maybe_string_can_pass_through_input(test_input):
    assert rn.maybe_str(test_input) == test_input


@pytest.mark.parametrize('test_input', ['', None])
def test_maybe_string_deals_with_missing_correctly(test_input):
    assert rn.maybe_str(test_input) is None


def test_maybe_string_can_strip_spaces():
    assert rn.maybe_str('  a  ', strip=True) == 'a'
    assert rn.maybe_str(1, strip=True) == '1'


@pytest.mark.parametrize('test_input', ['   ', None])
def test_maybe_string_can_strip_and_handle_missing(test_input):
    assert rn.maybe_str(test_input, strip=True) is None


@pytest.mark.parametrize('test_input', [1, '1', '1.0'])
def test_maybe_float_parses_floats(test_input):
    assert rn.maybe_float(test_input) == 1.0


@pytest.mark.parametrize('test_input', ['', None])
def test_maybe_float_handles_missing(test_input):
    assert rn.maybe_float(test_input) is None


def test_maybe_float_can_strip_whitespace():
    assert rn.maybe_float('   1   ') == 1
    assert rn.maybe_float('      ') is None


def test_maybe_int_parses_ints():
    assert rn.maybe_int(None) is None
    assert rn.maybe_int('') is None
    assert rn.maybe_int('1') == 1
    assert rn.maybe_int('   1    ', strip=True) == 1


@pytest.mark.parametrize('test_input', ['1', 1, '1    '])
def test_rename_build_useful_function(test_input):
    fn = rn.rename('a', rn.maybe_int, strip=True)
    assert fn({'a': test_input}) == 1
    assert fn.initial == 'a'


@pytest.mark.parametrize('test_input', ['', None, '    '])
def test_rename_handles_missing_or_none(test_input):
    fn = rn.rename('a', rn.maybe_int, strip=True)
    assert fn({'a': test_input}) is None
    assert fn.initial == 'a'
    assert fn({}) is None


def test_transform_keeps_names_same():
    fn = rn.transform('A-B', rn.maybe_int)
    assert fn({'A-B': '1', 'b': None}) == 1
    assert fn.final == 'A-B'


def test_with_dashes_cleans_up_names():
    fn = rn.with_dashes('a-B', rn.maybe_float, strip=True)
    assert fn({'a-B': '2.0  ', 'b': 1}) == 2.0
    assert fn.final == 'a_b'


class RenamerTest(ut.TestCase):
    def setUp(self):
        pass

    def test_can_use_renamers(self):
        renamer = rn.Renamer(pdb_id=rn.rename('pdb', str))
        assert renamer({'pdb': 'a', 'bob': 1}) == {'pdb_id': 'a'}
