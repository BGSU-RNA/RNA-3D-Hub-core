from cStringIO import StringIO
from datetime import datetime

from pymotifs import core
from pymotifs import utils
from pymotifs import models as mod


class Parser(object):
    def __call__(self, text):
        data = []
        for line in StringIO(text).readlines():
            # OBSLTE    26-SEP-06 2H33     2JM5 2OWI
            if 'OBSLTE' in line.split():
                parts = line.split()
                obsolete_date = datetime.strptime(parts[1], '%d-%b-%y')
                replaced = ",".join(parts[3:])
                if not replaced:
                    replaced = None
                data.append({
                    'pdb_obsolete_id': parts[2],
                    'date': obsolete_date,
                    'replaced_by': replaced
                })
        return data


class Loader(core.MassLoader):
    merge_data = True
    table = mod.PdbObsolete

    def has_data(self, pdb_id, **kwargs):
        return False

    def data(self, *args, **kwargs):
        """
        Download the file with all obsolete structures over ftp, store
        the data in the database, remove obsolete entries
        from the pdb_info table
        """

        try:
            ftp = utils.FTPFetchHelper('ftp.wwpdb.org', parser=Parser())
            return ftp('/pub/pdb/data/status/obsolete.dat')
        except Exception as err:
            self.logger.critical("Could not get obsolete ids")
            self.logger.exception(err)
            raise core.StageFailed("Could not get obsolete ids")
