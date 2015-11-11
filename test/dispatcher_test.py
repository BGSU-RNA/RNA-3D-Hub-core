import unittest as ut

from pymotifs.dispatcher import Dispatcher
from pymotifs.units import loader as units
from pymotifs.pdbs import loader as pdbs
from pymotifs import download as down
from pymotifs.export import loader as export
from pymotifs import mat_files as mat

from pymotifs.ife.info import Loader as IfeLoader
from pymotifs.interactions.pairwise import Loader as InteractionLoader


class DispatcherTest(ut.TestCase):
    def setUp(self):
        self.dispatcher = Dispatcher('units.info')

    def test_it_can_load_a_stage_by_name(self):
        val = self.dispatcher.get_stage('units.info')
        self.assertEquals(val, units.InfoLoader)

    def test_it_can_load_an_aggregate_by_name(self):
        val = self.dispatcher.get_stage('units.loader')
        self.assertEquals(val, units.Loader)

    def test_it_can_load_a_unit_stage(self):
        val = self.dispatcher.get_stage('units.quality')
        self.assertEquals(val, units.QualityLoader)

    def test_it_fails_if_given_unknown_module(self):
        self.assertRaises(ImportError, self.dispatcher.get_stage, 'bob')


class LoadingTest(ut.TestCase):
    def setUp(self):
        self.dispatcher = Dispatcher('')

    def test_it_knows_if_something_is_a_loader(self):
        loader = self.dispatcher.is_loader('pymotifs.units.info')
        val = loader(units.InfoLoader)
        self.assertTrue(val)

    def test_it_can_load_a_single_loader_from_many_imported(self):
        loader = self.dispatcher.is_loader('pymotifs.units.loader')
        val = loader(units.Loader)
        self.assertTrue(val)

    def test_it_will_skip_those_with_mismatched_name(self):
        loader = self.dispatcher.is_loader('pymotifs.units.info')
        val = loader(units.Loader)
        self.assertFalse(val)


class LoadingStagesTest(ut.TestCase):
    def setUp(self):
        self.dispatcher = Dispatcher('units.info')

    def test_can_give_only_selected_stage(self):
        self.dispatcher.skip_dependencies = True
        val = self.dispatcher.stages('units.info')
        self.assertEquals([units.InfoLoader], val)

    def test_can_get_with_all_dependecies(self):
        val = self.dispatcher.stages('units.info')
        ans = [down.Downloader, pdbs.InfoLoader, units.InfoLoader]
        self.assertEquals(ans, val)

    def test_can_get_all_dependencies_of_multistage(self):
        val = self.dispatcher.stages('units.loader')
        ans = [down.Downloader, pdbs.InfoLoader,
               units.InfoLoader, export.CifAtom,
               units.DistancesLoader, units.QualityLoader,
               mat.Loader, units.RedundantNucleotidesLoader]
        self.assertEquals(ans, val)

    def test_can_exclude_specific_stages(self):
        self.dispatcher.exclude = set(['units.distances'])
        val = self.dispatcher.stages('units.loader')
        ans = [down.Downloader, pdbs.InfoLoader,
               units.InfoLoader, export.CifAtom,
               units.QualityLoader, mat.Loader,
               units.RedundantNucleotidesLoader]
        self.assertEquals(ans, val)

    def test_it_respects_ordering(self):
        stages = self.dispatcher.stages('update')
        index1 = stages.index(InteractionLoader)
        index2 = stages.index(IfeLoader)
        self.assertTrue(index1 < index2)

    def test_can_exclude_a_stage_collection(self):
        self.dispatcher.exclude = set(['pdbs.loader', 'units.distances'])
        val = self.dispatcher.stages('units.loader')
        ans = [down.Downloader,
               units.InfoLoader, export.CifAtom,
               units.QualityLoader, mat.Loader,
               units.RedundantNucleotidesLoader]
        self.assertEquals(ans, val)
