from pymotifs import core
from pymotifs import models as mod
from pymotifs.ss.helpers import Parser
from pymotifs.ss.helpers import md5

from pymotifs.ss.info import Loader as SsInfoLoader


class Loader(core.SimpleLoader):
    dependencies = set([SsInfoLoader])
    @property
    def table(self):
        return mod.SsPositions
    mark = False

    def ss_id(self, filename):
        with self.session() as session:
            return session.query(mod.SsInfo).\
                filter_by(md5=md5(filename)).\
                one().\
                ss_id

    def to_process(self, pdbs, **kwargs):
        if 'filename' not in kwargs:
            raise core.Skip("No filename to import")
        return [kwargs['filename']]

    def query(self, session, filename):
        return session.query(self.table).\
            join(mod.SsInfo, mod.SsInfo.ss_id == self.table.ss_id).\
            filter(mod.SsInfo.md5 == md5(filename))

    def data(self, ss_filename, **kwargs):
        parser = Parser(self.config, self.session)
        parsed = parser(ss_filename)
        nts = parsed['nts']
        ss_id = self.ss_id(ss_filename)
        for nt in nts:
            nt['ss_id'] = ss_id
        return nts
