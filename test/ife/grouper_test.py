from test import StageTest

from pymotifs.ife.grouper import Grouper as AutonomousGrouper


class AutonomousCheckTest(StageTest):
    loader_class = AutonomousGrouper

    def test_requires_majority_internal_bp(self):
        chain = {'internal': 30, 'external': 30}
        self.assertTrue(self.loader.is_autonomous(chain))

    def test_is_not_autonomous_if_few_interactions(self):
        chain = {'internal': 3, 'external': 3}
        self.assertFalse(self.loader.is_autonomous(chain))

    def test_knows_if_not_autonomous(self):
        chain = {'internal': 30, 'external': 100}
        self.assertFalse(self.loader.is_autonomous(chain))

    def test_allows_no_internal_bp(self):
        chain = {'internal': 0, 'external': 10}
        self.assertFalse(self.loader.is_autonomous(chain))

    def test_allows_no_external_bps(self):
        chain = {'internal': 10, 'external': 0}
        self.assertTrue(self.loader.is_autonomous(chain))


class GroupingTest(StageTest):
    loader_class = AutonomousGrouper

    def test_will_mark_chains_as_autonomous(self):
        chains = [{'internal': 10, 'external': 0, 'pdb': '1S72', 'name': 'A',
                   'id': '1S72|1|A'}]
        ans = [[{'internal': 10, 'external': 0, 'pdb': '1S72', 'name': 'A',
                'id': '1S72|1|A', 'autonomous': True}]]
        val = self.loader.group(chains, {})
        self.assertEquals(ans, val)

    def test_will_put_non_interacting_chains_seperately(self):
        chains = [{'internal': 10, 'external': 0, 'pdb': '1S72', 'name': 'A',
                   'id': '1S72|1|A'},
                  {'internal': 10, 'external': 0, 'pdb': '1S72', 'name': 'B',
                   'id': '1S71|1|B'}]
        ans = [[{'internal': 10, 'external': 0, 'pdb': '1S72', 'name': 'A',
                'id': '1S72|1|A', 'autonomous': True}],
               [{'internal': 10, 'external': 0, 'pdb': '1S72', 'name': 'B',
                'id': '1S71|1|B', 'autonomous': True}]]
        val = self.loader.group(chains, {})
        self.assertEquals(ans, val)

    def test_puts_non_interactions_non_autonomous_chains_seperately(self):
        chains = [{'internal': 0, 'external': 10, 'pdb': '1S72', 'name': 'A',
                   'id': '1S72|1|A'},
                  {'internal': 0, 'external': 10, 'pdb': '1S72', 'name': 'B',
                   'id': '1S71|1|B'}]
        ans = [[{'internal': 0, 'external': 10, 'pdb': '1S72', 'name': 'A',
                'id': '1S72|1|A', 'autonomous': False}],
               [{'internal': 0, 'external': 10, 'pdb': '1S72', 'name': 'B',
                'id': '1S71|1|B', 'autonomous': False}]]
        val = self.loader.group(chains, {})
        self.assertEquals(ans, val)

    def test_will_put_autonomous_but_interacting_chains_seperatly(self):
        chains = [{'internal': 10, 'external': 0, 'pdb': '1S72', 'name': 'A',
                   'id': '1S72|1|A'},
                  {'internal': 10, 'external': 0, 'pdb': '1S72', 'name': 'B',
                   'id': '1S71|1|B'}]
        ans = [[{'internal': 10, 'external': 0, 'pdb': '1S72', 'name': 'A',
                'id': '1S72|1|A', 'autonomous': True}],
               [{'internal': 10, 'external': 0, 'pdb': '1S72', 'name': 'B',
                'id': '1S71|1|B', 'autonomous': True}]]
        interactions = {'A': {'B': 10}, 'B': {'A': 10}}
        val = self.loader.group(chains, interactions)
        self.assertEquals(ans, val)

    def test_will_join_non_autonomous_chains_on_interactions(self):
        chains = [{'internal': 0, 'external': 10, 'pdb': '1S72', 'name': 'A',
                   'id': '1S72|1|A'},
                  {'internal': 0, 'external': 10, 'pdb': '1S72', 'name': 'B',
                   'id': '1S71|1|B'}]
        ans = [[{'internal': 0, 'external': 10, 'pdb': '1S72', 'name': 'A',
                 'id': '1S72|1|A', 'autonomous': False},
                {'internal': 0, 'external': 10, 'pdb': '1S72', 'name': 'B',
                 'id': '1S71|1|B', 'autonomous': False}]]
        interactions = {'A': {'B': 10}, 'B': {'A': 10}}
        val = self.loader.group(chains, interactions)
        self.assertEquals(ans, val)

    def test_non_automous_single_chains_will_be_alone(self):
        chains = [{'internal': 0, 'external': 10, 'pdb': '1S72', 'name': 'A',
                   'id': '1S72|1|A'}]
        ans = [[{'internal': 0, 'external': 10, 'pdb': '1S72', 'name': 'A',
                'id': '1S72|1|A', 'autonomous': False}]]
        val = self.loader.group(chains, {})
        self.assertEquals(ans, val)


class SplittingGroupsByAutonomous(StageTest):
    """The idea here is that one or more non-autonomous chains will form a
    group with more than one autonomous chain. In this case we want the
    autonomous chains to be alone along with whatever should accompany them.
    These tests try to confirm that.
    """
    loader_class = AutonomousGrouper

    def test_will_split_two_into_two_groups(self):
        chains = [{'internal': 10, 'external': 0, 'pdb': '1S72', 'name': 'A',
                   'id': '1S72|1|A', 'autonomous': True},
                  {'internal': 10, 'external': 0, 'pdb': '1S72', 'name': 'B',
                   'id': '1S71|1|B', 'autonomous': True},
                  {'internal': 0, 'external': 10, 'pdb': '1S72', 'name': 'C',
                   'id': '1S71|1|C', 'autonomous': False}]
        ans = [
            [{'internal': 10, 'external': 0, 'pdb': '1S72', 'name': 'A',
              'id': '1S72|1|A', 'autonomous': True},
             {'internal': 0, 'external': 10, 'pdb': '1S72', 'name': 'C',
              'id': '1S71|1|C', 'autonomous': False}],
            [{'internal': 10, 'external': 0, 'pdb': '1S72', 'name': 'B',
              'id': '1S71|1|B', 'autonomous': True},
             {'internal': 0, 'external': 10, 'pdb': '1S72', 'name': 'C',
              'id': '1S71|1|C', 'autonomous': False}]]
        interactions = {'A': {'C': 10},
                        'B': {'C': 10},
                        'C': {'A': 10, 'B': 10}}
        val = self.loader.split_group(chains, interactions)
        self.assertEquals(ans, val)

    def test_grouping_includes_group_splits(self):
        chains = [{'internal': 10, 'external': 0, 'pdb': '1S72', 'name': 'A',
                   'id': '1S72|1|A', 'autonomous': True},
                  {'internal': 10, 'external': 0, 'pdb': '1S72', 'name': 'B',
                   'id': '1S71|1|B', 'autonomous': True},
                  {'internal': 0, 'external': 10, 'pdb': '1S72', 'name': 'C',
                   'id': '1S71|1|C', 'autonomous': False}]
        ans = [
            [{'internal': 10, 'external': 0, 'pdb': '1S72', 'name': 'A',
              'id': '1S72|1|A', 'autonomous': True},
             {'internal': 0, 'external': 10, 'pdb': '1S72', 'name': 'C',
              'id': '1S71|1|C', 'autonomous': False}],
            [{'internal': 10, 'external': 0, 'pdb': '1S72', 'name': 'B',
              'id': '1S71|1|B', 'autonomous': True},
             {'internal': 0, 'external': 10, 'pdb': '1S72', 'name': 'C',
              'id': '1S71|1|C', 'autonomous': False}]]
        interactions = {'A': {'C': 10},
                        'B': {'C': 10},
                        'C': {'A': 10, 'B': 10}}
        val = self.loader.group(chains, interactions)
        self.assertEquals(ans, val)


class UsingCorrectChainsTest(StageTest):
    loader_class = AutonomousGrouper

    def test_uses_correct_chains(self):
        chains = self.loader('1FEU')
        val = [c['id'] for c in chains]
        ans = ['1FEU|C+1FEU|B', '1FEU|F+1FEU|E']
        self.assertEquals(ans, val)
