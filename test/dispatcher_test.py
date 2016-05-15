import pytest

import unittest as ut

from test import CONFIG
from test import Session

from pymotifs import core
from pymotifs.dispatcher import Dispatcher
from pymotifs.units import loader as units
from pymotifs.pdbs import loader as pdbs
from pymotifs.ife import loader as ifes
from pymotifs import download as down
from pymotifs.chains import loader as chains
from pymotifs.interactions import loader as interactions
# from pymotifs.export import loader as export
# from pymotifs import mat_files as mat
from pymotifs.cli.introspect import UnknownStageError

# from pymotifs.ife.info import Loader as IfeLoader
# from pymotifs.interactions.pairwise import Loader as InteractionLoader


class ToExcludeTests(ut.TestCase):
    def setUp(self):
        self.dispatcher = Dispatcher('units.info', CONFIG, Session)

    def test_given_nothing_returns_nothing(self):
        assert self.dispatcher.to_exclude() == set()

    def test_given_skip_dependencies_returns_true(self):
        self.dispatcher.skip_dependencies = True
        assert self.dispatcher.to_exclude(skip_dependencies=True) is True
        val = self.dispatcher.to_exclude('units.info', skip_dependencies=True)
        assert val is True

    def test_will_dedup_excludes(self):
        val = self.dispatcher.to_exclude('units.info', 'units.info')
        assert val == set(['units.info'])

    def test_will_complain_given_unknown_stage(self):
        with pytest.raises(UnknownStageError):
            self.dispatcher.to_exclude('bob')

    def test_will_expand_stage_container_to_stages(self):
        val = self.dispatcher.to_exclude('units.loader')
        assert val == set(['units.info', 'units.distances', 'units.quality',
                           'units.distances', 'units.centers', 'units.loader',
                           'units.rotation', 'units.redundant'])


class DependenciesTest(ut.TestCase):
    def setUp(self):
        self.dispatcher = Dispatcher('units.info', CONFIG, Session)

    def test_it_can_compute_a_dependency_graph(self):
        val = self.dispatcher.dependencies([units.InfoLoader])
        assert val == {
            units.InfoLoader: set([down.Downloader, pdbs.InfoLoader]),
            pdbs.InfoLoader: set([]),
            down.Downloader: set([]),
        }

    def test_it_computes_dependencies_for_stage_container(self):
        val = self.dispatcher.dependencies([pdbs.Loader])
        assert val == {
            pdbs.InfoLoader: set([]),
            pdbs.ObsoleteLoader: set([pdbs.InfoLoader]),
        }

    def test_it_computes_dependencies_for_several(self):
        val = self.dispatcher.dependencies([pdbs.Loader, units.InfoLoader])
        assert val == {
            pdbs.InfoLoader: set([]),
            down.Downloader: set([]),
            pdbs.ObsoleteLoader: set([pdbs.InfoLoader]),
            units.InfoLoader: set([down.Downloader, pdbs.InfoLoader]),
        }

    def test_it_computes_dependencies_for_nested_containers(self):
        val = self.dispatcher.dependencies([ifes.Loader])
        assert val[ifes.InfoLoader] == \
            set([chains.InfoLoader, chains.SpeciesLoader,
                 interactions.PairwiseLoader, interactions.HLSummaryLoader])


class LevelsTest(ut.TestCase):
    def setUp(self):
        self.dispatcher = Dispatcher('units.info', CONFIG, Session)

    def levels(self, *args):
        levels = []
        for objs in self.dispatcher.levels(*args):
            levels.append([obj.name for obj in objs])
        return levels

    def test_produces_correct_number_of_levels(self):
        deps = self.dispatcher.dependencies([pdbs.Loader])
        val = self.levels(deps, [], True)
        assert val == [['pdbs.info'], ['pdbs.obsolete']]

    def test_has_no_excluded_stages(self):
        deps = self.dispatcher.dependencies([pdbs.Loader])
        val = self.levels(deps, ['pdbs.obsolete'], [])
        assert ['pdbs.obsolete'] not in val
        assert val == [['pdbs.info']]

    def test_will_always_have_allowed_stages(self):
        deps = self.dispatcher.dependencies([pdbs.Loader])
        val = self.levels(deps, ['pdbs.obsolete'], ['pdbs.obsolete'])
        assert val == [['pdbs.info'], ['pdbs.obsolete']]


class StagesTest(ut.TestCase):
    def setUp(self):
        self.dispatcher = Dispatcher('units.info', CONFIG, Session)

    def stages(self, *args):
        return [o.name for o in self.dispatcher.stages(*args)]

    def test_can_load_requested_stages(self):
        assert self.stages('units.info') == [
            'download',
            'pdbs.info',
            'units.info',
        ]

    def test_when_skipping_produces_container_stages_in_correct_order(self):
        self.dispatcher.skip_dependencies = True
        assert self.stages('units.loader') == [
            'units.info',
            'units.centers',
            'units.distances',
            'units.quality',
            'units.rotation',
            'units.redundant',
        ]

    def test_it_can_load_stage_container(self):
        assert self.stages('units.loader') == [
            'download',
            'pdbs.info',
            'export.cifatom',
            'units.info',
            'mat_files',
            'units.centers',
            'units.distances',
            'units.quality',
            'units.rotation',
            'units.redundant',
        ]

    def test_can_exclude_specific_stages(self):
        self.dispatcher.exclude = set(['units.distances'])
        val = self.stages('units.loader')
        assert 'units.distances' not in val
        assert val == [
            'download',
            'pdbs.info',
            'export.cifatom',
            'units.info',
            'mat_files',
            'units.centers',
            'units.quality',
            'units.rotation',
            'units.redundant',
        ]

    def test_can_exclude_a_stage_collection(self):
        self.dispatcher.exclude = set(['pdbs.loader', 'units.distances'])
        assert self.stages('units.loader') == [
            'download',
            'export.cifatom',
            'units.info',
            'mat_files',
            'units.centers',
            'units.quality',
            'units.rotation',
            'units.redundant',
        ]

    @pytest.mark.xfail(reason="Haven't worked on yet")
    def test_complains_if_nothing_to_run(self):
        self.dispatcher.skip_dependencies = True
        self.dispatcher.exclude.add('units.info')
        with pytest.raises(core.InvalidState):
            self.stages('units.info')

    def test_can_get_one_stage_built(self):
        self.dispatcher.skip_dependencies = True
        assert self.stages('units.info') == ['units.info']

    def test_running_with_nested_container_dependencies(self):
        val = self.stages('ife.loader')
        assert val.index('interactions.pairwise') < val.index('ife.info')
        assert val == [
            'download',
            'pdbs.info',
            'export.cifatom',
            'pdbs.obsolete',
            'units.info',
            'chains.info',
            'mat_files',
            'interactions.pairwise',
            'species_mapping',
            'chains.species',
            'interactions.helix_loop_summary',
            'ife.info'
        ]
