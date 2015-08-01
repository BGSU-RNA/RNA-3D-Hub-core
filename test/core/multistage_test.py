from test import StageTest

from pymotifs import core


class C(core.Stage):
    pass


class D(core.Stage):
    pass


class A(core.Stage):
    dependencies = set([C])


class B(core.Stage):
    dependencies = set([C, D])


class Multi(core.MultiStageLoader):
    stages = [A, B]


class DependenciesTest(StageTest):
    loader_class = Multi

    def test_it_computes_all_dependencies_of_stage_instances(self):
        self.assertEqual(set([C, D]), self.loader.dependencies)

    def test_it_computes_all_dependencies_of_stage_class(self):
        self.assertEqual(set([C, D]), self.loader_class.dependencies)
