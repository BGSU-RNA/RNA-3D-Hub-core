import unittest as ut

from pymotifs.nr.connectedsets import findconnectedsets as conn


class ConnectionTest(ut.TestCase):

    def test_returns_empty_given_none(self):
        self.assertEquals({}, conn({}))
