import unittest as ut

from test import CONFIG
from test import Session

from pymotifs.dispatcher import Dispatcher
from pymotifs.units import loader as units
from pymotifs.pdbs import loader as pdbs
from pymotifs import download as down
from pymotifs.export import loader as export
from pymotifs import mat_files as mat

from pymotifs.ife.info import Loader as IfeLoader
from pymotifs.interactions.pairwise import Loader as InteractionLoader


class LoadingStagesTest(ut.TestCase):
    def setUp(self):
        self.dispatcher = Dispatcher('units.info', CONFIG, Session)

    def test_can_give_only_selected_stage(self):
        self.dispatcher.skip_dependencies = True
        val = self.dispatcher.stages('units.info')
        self.assertEquals([units.InfoLoader], val)

    def test_given_mass_gets_all_stages_when_skipping_dependencies(self):
        self.dispatcher.skip_dependencies = True
        val = self.dispatcher.stages('units.loader')
        ans = [units.InfoLoader, units.QualityLoader, units.DistancesLoader,
               units.RedundantNucleotidesLoader]
        self.assertEquals(ans, val)

    def test_can_get_one_stage_built(self):
        self.dispatcher.skip_dependencies = True
        val = self.dispatcher.stages('units.info', build=True)
        self.assertEquals(1, len(val))
        self.assertTrue(isinstance(val[0], units.InfoLoader))

    def test_can_get_with_all_dependecies(self):
        val = self.dispatcher.stages('units.info')
        ans = [down.Downloader, pdbs.InfoLoader, units.InfoLoader]
        self.assertEquals(ans, val)

    def test_can_get_all_dependencies_of_multistage(self):
        val = self.dispatcher.stages('units.loader')
        ans = [down.Downloader, pdbs.InfoLoader,
               export.CifAtom, units.InfoLoader,
               mat.Loader, units.DistancesLoader,
               units.QualityLoader, units.RedundantNucleotidesLoader]
        self.assertEquals(ans, val)

    def test_can_exclude_specific_stages(self):
        self.dispatcher.exclude = set(['units.distances'])
        val = self.dispatcher.stages('units.loader')
        ans = [down.Downloader, pdbs.InfoLoader,
               export.CifAtom, units.InfoLoader,
               mat.Loader, units.QualityLoader,
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
        ans = [down.Downloader, export.CifAtom,
               units.InfoLoader,
               mat.Loader, units.QualityLoader,
               units.RedundantNucleotidesLoader]
        self.assertEquals(ans, val)
