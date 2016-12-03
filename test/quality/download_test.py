import os

from test import StageTest

from pymotifs.quality.download import Loader

import pytest


class DownloaderTest(StageTest):
    loader_class = Loader

    def test_it_can_detect_if_has_no_data(self):
        assert self.loader.has_data('0FJG') is False

    def test_knows_if_has_data(self):
        if not self.loader.has_data('1FJG'):
            self.loader(['1FJG'])
        assert self.loader.has_data('1FJG') is True

    @pytest.mark.skip()
    def test_can_remove_data(self):
        pass

    def test_can_get_filename(self):
        ans = os.path.abspath('MotifAtlas/quality/validation-reports/1FJG.xml.gz')
        assert self.loader.filename('1FJG') == ans

    def test_gives_empty_string_for_missing(self):
        assert self.loader.data('0FJG') == ''

    def test_can_fetch_data(self):
        val = self.loader.data('1S72')
        assert len(val) > 0

    def test_can_download_a_file(self):
        if not self.loader.has_data('1FJG'):
            self.loader.remove('1FJG')

        assert os.path.exists(self.loader.filename('1FJG')) is False
        self.loader(['1FJG'])
        assert os.path.exists(self.loader.filename('1FJG')) is True
