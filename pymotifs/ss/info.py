import hashlib

from pymotifs import core
from pymotifs import models as mod
from pymotifs.ss.helpers import Parser


class Loader(core.SimpleLoader):
    table = mod.SsInfo
    mark = False

    def md5(self, filename):
        with open(filename, 'rb') as raw:
            return hashlib.md5(raw.read).hexdigest()

    def to_process(self, pdbs, **kwargs):
        if 'filename' not in kwargs:
            raise core.Skip("No filename to import")
        return [kwargs['filename']]

    def query(self, session, filename):
        return session.query(self.table).filter_by(md5=self.md5(filename))

    def data(self, filename, **kwargs):
        parsed = Parser(filename)
        return {
            'name': parsed['name'],
            'filename': filename,
            'sequence': parsed['sequence'],
            'md5': self.md5(filename)
        }
