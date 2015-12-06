from pymotifs import core
from pymotifs import models as mod
from pymotifs.ss.helpers import Parser
from pymotifs.ss.helpers import md5


class Loader(core.SimpleLoader):
    table = mod.SsInfo
    mark = False

    def to_process(self, pdbs, **kwargs):
        if 'filename' not in kwargs:
            raise core.Skip("No filename to import")
        return [kwargs['filename']]

    def query(self, session, filename):
        return session.query(self.table).filter_by(md5=md5(filename))

    def data(self, ss_filename, **kwargs):
        parser = Parser(self.config, self.session)
        parsed = parser(ss_filename)
        return {
            'name': kwargs.get('ss_name', None),
            'filename': ss_filename,
            'sequence': parsed['sequence'],
            'md5': md5(ss_filename)
        }
