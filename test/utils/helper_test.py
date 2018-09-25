import unittest

import pymotifs.utils as utils


class FailingRetry(utils.RetryHelper):
    def action(*args, **kwargs):
        raise ValueError("Expected")


class RetryHelperTest(unittest.TestCase):
    def test_can_allow_failure(self):
        val = FailingRetry(allow_fail=True)()
        self.assertTrue(val is None)

    def test_will_fail_with_failure(self):
        self.assertRaises(utils.RetryFailedException, FailingRetry())


class GzipFetchHelperTest(unittest.TestCase):
    def setUp(self):
        self.helper = utils.GzipFetchHelper()

    def test_can_fetch_a_zipped(self):
        val = self.helper('http://www.rcsb.org/pdb/files/1GID.pdb.gz')
        self.assertTrue(val is not None)

    def test_unzips_fetched(self):
        val = self.helper('http://www.rcsb.org/pdb/files/1GID.pdb.gz')
        self.assertEquals('HEADER', val[0:6])
