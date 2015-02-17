import os

from test import StageTest
from test import skip_without_matlab
from nose import SkipTest

from pymotifs.loops.quality import LoopQualityLoader
from pymotifs.models import LoopQa


class ParsingTest(StageTest):
    loader_class = LoopQualityLoader

    def setUp(self):
        super(ParsingTest, self).setUp()
        with open('test/files/loop-quality.csv', 'rb') as raw:
            self.data = self.loader.parse(raw, '1.9')

    def test_can_parse_all_of_a_file(self):
        self.assertEquals(20, len(self.data))

    def test_creates_correct_data(self):
        ans = {
            'id': "HL_1GID_001",
            'status': "1",
            'modifications': None,
            'nt_signature': "149, 150, 151, 152, 153, 154",
            'complementary': None,
            'release_id': '1.9'
        }
        self.assertEquals(ans, self.data[0])


class QueryingTest(StageTest):
    loader_class = LoopQualityLoader

    def test_can_detect_if_has_data(self):
        self.assertTrue(self.loader.has_data(('1GID', '1.0')))

    def test_can_detect_if_has_no_data(self):
        self.assertFalse(self.loader.has_data(('0GID', '1.0')))

    def test_can_remove_known_data(self):
        raise SkipTest()

    def test_works_if_removing_missing_data(self):
        raise SkipTest()


@skip_without_matlab
class MatlabTests(StageTest):
    loader_class = LoopQualityLoader

    def setUp(self):
        super(ParsingTest, self).setUp()
        self.data = self.loader.data(('1GID', '2.0'))

    def test_gets_all_loop_quality(self):
        self.assertEquals(22, len(self.data))

    def test_creates_correct_data(self):
        ans = LoopQa(id="HL_1GID_001",
                     status="1",
                     modifications=None,
                     nt_signature="149, 150, 151, 152, 153, 154",
                     complementary=None,
                     release_id='2.0')
        self.assertEquals(ans, self.data[0])

    def test_removes_data_file(self):
        self.assertFalse(os.path.exists('./FR3D/LoopQA.csv'))
