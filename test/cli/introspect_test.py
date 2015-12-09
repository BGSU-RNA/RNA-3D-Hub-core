import test  # Used to make sure all models are reflected
from unittest import TestCase

from pymotifs.cli import introspect as intro

from pymotifs.ife import info
from pymotifs.units import loader as units
from pymotifs import update


class GettingStageTest(TestCase):
    def test_can_get_a_stage_by_name(self):
        self.assertEquals(info, intro.get_stage('ife.info'))

    def test_complains_given_unknown_name(self):
        self.assertRaises(intro.UnknownStageError, intro.get_stage, 'bob')

    def test_can_get_documentation_for_a_stage(self):
        val = intro.get_stage_info('ife.info')
        self.assertEquals('ife.info', val[0])
        self.assertEquals('Load IFE data', val[1])
        self.assertTrue(val[2])

    def test_can_determine_if_stage_has_loaders(self):
        self.assertTrue(intro.has_stage_loader(info))

    def test_can_determine_if_stage_has_no_loaders(self):
        self.assertFalse(intro.has_stage_loader(test))

    def test_can_determine_if_is_stage(self):
        self.assertTrue(intro.is_stage('ife.info'))

    def test_can_determine_if_is_not_stage(self):
        self.assertFalse(intro.is_stage('ife.helpers'))

    def test_it_can_load_a_stage_by_name(self):
        val = intro.get_stage('units.info')
        self.assertEquals(val, units.InfoLoader)

    def test_it_can_load_an_aggregate_by_name(self):
        val = intro.get_stage('units.loader')
        self.assertEquals(val, units.Loader)

    def test_it_can_load_a_unit_stage(self):
        val = intro.get_stage('units.quality')
        self.assertEquals(val, units.QualityLoader)


class GettingLoadersTest(TestCase):
    def test_can_get_a_loader_by_stage_name(self):
        self.assertEquals(info.Loader, intro.get_loader('ife.info'))

    def test_complains_given_unknown_name(self):
        self.assertRaises(intro.UnknownStageError, intro.get_stage, 'bob')

    def test_can_detect_if_is_stage_loader(self):
        self.assertTrue(intro.is_stage_loader(info, info.Loader))

    def test_can_detect_if_is_not_stage_loader(self):
        self.assertFalse(intro.is_stage_loader(test, info.Loader))

    def test_can_get_a_mass_loader(self):
        self.assertEquals(update.Loader, intro.get_loader('update'))
