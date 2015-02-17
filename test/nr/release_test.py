from test import StageTest

from pymotifs.nr.release import Loader


class ReleaseTest(StageTest):
    loader_class = Loader

    def test_gets_next_release_id(self):
        self.loader.config['release_mode']['nrlist'] = 'minor'
        data = self.loader.data()
        self.assertEquals('1.71', data.id)

    def test_gets_next_release_using_config(self):
        self.loader.config['release_mode']['nrlist'] = 'major'
        data = self.loader.data()
        self.assertEquals('2.0', data.id)
