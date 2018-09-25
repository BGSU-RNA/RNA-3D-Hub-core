from test import StageTest

from pymotifs.nr.release import Loader


class ReleaseTest(StageTest):
    loader_class = Loader

    def test_gets_next_release_id(self):
        self.loader.config['release_mode']['nrlist'] = 'minor'
        assert self.loader.next_id('3.0') == '3.1'

    def test_gets_next_release_using_config(self):
        self.loader.config['release_mode']['nrlist'] = 'major'
        assert self.loader.next_id('3.0') == '4.0'
