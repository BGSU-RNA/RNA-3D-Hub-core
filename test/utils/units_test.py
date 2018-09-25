from unittest import TestCase

from fr3d.data import Component

from pymotifs.utils import units


class ComputingTheType(TestCase):

    def test_can_determine_is_rna(self):
        component = Component([], sequence='A')
        self.assertEquals('rna', units.component_type(component))

    def test_can_determine_is_dna(self):
        component = Component([], sequence='DA')
        self.assertEquals('dna', units.component_type(component))

    def test_can_determine_is_aa(self):
        component = Component([], sequence='arg')
        self.assertEquals('aa', units.component_type(component))

    def test_can_determine_is_water(self):
        component = Component([], sequence='hoh')
        self.assertEquals('water', units.component_type(component))

    def test_gives_none_otherwise(self):
        component = Component([], sequence='gtp')
        self.assertEquals(None, units.component_type(component))
