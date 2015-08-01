import unittest as ut

from pymotifs.utils.connectedsets import find_connected as conn


class ConnectionTest(ut.TestCase):

    def test_returns_empty_given_none(self):
        self.assertEquals({}, conn({}))

    def test_basic_test(self):
        connections = {}
        connections['A'] = ['B', 'C']
        connections['B'] = ['D']
        connections['C'] = ['E']
        connections['E'] = ['A']
        connections['F'] = ['B']
        connections['zA'] = ['zB', 'zC']
        connections['zB'] = ['zD']
        connections['zC'] = ['zE']
        connections['zD'] = ['zF']
        connections['zE'] = ['zA']
        connections['zF'] = ['zB']
        for key in connections.keys():
            connections[key] = set(connections[key])
        ans = {'A': set(['A', 'C', 'B', 'E', 'D', 'F']),
               'zD': set(['zD', 'zE', 'zF', 'zA', 'zB', 'zC'])}
        val = conn(connections)
        self.assertEquals(ans, val)
