from cStringIO import StringIO
from datetime import datetime

from pymotifs import core
from pymotifs import utils
from pymotifs.models import PdbObsolete


class Parser(object):
    def __call__(self, text):
        data = []
        for line in StringIO(text).readlines():
            # OBSLTE    26-SEP-06 2H33     2JM5 2OWI
            if 'OBSLTE' in line.split():
                parts = line.split()
                obsolete_date = datetime.strptime(parts[1], '%d-%b-%y')
                data.append({
                    'id': parts[2],
                    'date': obsolete_date,
                    'replaced_by': parts[3:]
                })
        return data


class Loader(core.MassLoader):
    merge_data = True

    def data(self, *args, **kwargs):
        """
        Download the file with all obsolete structures over ftp, store
        the data in the database, remove obsolete entries
        from the pdb_info table
        """

        try:
            ftp = utils.FTPFetchHelper('ftp.wwpdb.org', parser=Parser())
            parsed = ftp('/pub/pdb/data/status/obsolete.dat')
        except Exception as err:
            self.logger.critical("Could not get obsolete ids")
            self.logger.exception(err)
            raise core.StageFailed("Could not get obsolete ids")

        data = []
        for entry in parsed:
            replaced = ','.join(entry['replaced_by'])
            data.append(PdbObsolete(obsolete_id=entry['id'],
                                    date=entry['date'],
                                    replaced_by=replaced))

        if not data:
            self.logger.error("Found no obsolete ids.")

        return data
