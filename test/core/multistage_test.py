from test import StageTest

from pymotifs import core

RAN = []


class C(core.Stage):
    pass


class D(core.Stage):
    pass


class A(core.Stage):
    dependencies = set([C])

    def __call__(self, data, **kwargs):
        RAN.append(('A', data))


class B(core.Stage):
    dependencies = set([A, C, D])

    def __call__(self, data, **kwargs):
        RAN.append(('B', data))


class Multi(core.MultiStageLoader):
    stages = set([A, B])


class DependenciesTest(StageTest):
    loader_class = Multi

    def setUp(self):
        super(DependenciesTest, self).setUp()
        RAN = []

    def test_it_computes_all_dependencies_of_stage_instances(self):
        self.assertEqual(set([C, D]), self.loader.dependencies)

    def test_it_computes_all_dependencies_of_stage_class(self):
        self.assertEqual(set([C, D]), self.loader_class.dependencies)

    def test_runs_stages_in_order(self):
        self.loader(['data'])
        self.assertEqual(RAN, [('A', ['data']), ('B', ['data'])])
