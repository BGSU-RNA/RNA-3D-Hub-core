from test import StageTest

from pymotifs.models import ExpSeqInfo as Exp
from pymotifs.exp_seq.positions import Loader


class GeneratingDataTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(GeneratingDataTest, self).setUp()
        self.data = self.loader.positions(10, 'ACGU')

    def test_can_find_all_positions(self):
        self.assertEquals(4, len(self.data))

    def test_can_generate_valid_data(self):
        ans = {
            'exp_seq_id': 10,
            'unit': 'A',
            'index': 0
        }
        self.assertEquals(ans, self.data[0])


class GettingSequencesTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(GettingSequencesTest, self).setUp()
        self.loader.store([Exp(id=-1, sequence='ACGU', length=3),
                           Exp(id=-2, sequence='AAAAAAA', length=6)])
        self.data = sorted(self.loader.sequences())

    def tearDown(self):
        with self.loader.session() as session:
            session.query(Exp).\
                filter(Exp.id.in_([-1, -2])).\
                delete(synchronize_session=False)

    def test_finds_all_sequences(self):
        self.assertEquals(2, len(self.data))

    def test_gets_correct_data(self):
        self.assertEquals((-1, 'ACGU'), self.data[-1])

# class BuildingDataTest(Stag
