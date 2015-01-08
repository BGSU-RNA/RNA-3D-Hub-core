import core
import utils


class Loader(core.Stage):
    def __init__(self, *args):
        super(Loader, self).__init__(*args)
        self.ftp = utils.FTPFetchHelper('ftp.wwpdb.org')

    def __call__(self, *args, **kwargs):
        """Download the file with all obsolete structures over ftp, store
           the data in the database, remove obsolete entries
           from the pdb_info table"""

        TEMPFILE = 'obsolete.dat'
        REPEAT = 10
        done = False

        """download the data file from PDB"""
        for i in xrange(REPEAT):
            try:
                ftp = FTP('ftp.wwpdb.org')
                ftp.login()
                ftp.cwd('/pub/pdb/data/status')
                ftp.retrbinary("RETR %s" % TEMPFILE, open(TEMPFILE,"wb").write)
                ftp.quit()
                done = True
                logger.info('Downloaded obsolete.dat')
                break
            except Exception, e:
                logger.warning(e)
                logger.warning('Ftp download failed. Retrying...')

        if not done:
            logger.critical('All attempts to download obsolete.dat over ftp failed')
            logger.critical('Obsolete PDB files not updated')
            return

        """parse the data file"""
        obsolete_ids = []
        f = open(TEMPFILE, 'r')
        for line in f:
            # OBSLTE    26-SEP-06 2H33     2JM5 2OWI
            if 'OBSLTE' in line:
                parts = line.split()
                obsolete_date = datetime.strptime(parts[1], '%d-%b-%y')
                obsolete_ids.append(parts[2])
                replaced_by = ','.join(parts[3:])
                dbObj = PdbObsolete(obsolete_id=obsolete_ids[-1],
                                    date=obsolete_date,
                                    replaced_by=replaced_by)
                session.merge(dbObj)
        session.commit()

        """remove obsoleted files from pdb_info"""
        session.query(PdbInfo).\
                filter(PdbInfo.structureId.in_(obsolete_ids)).\
                delete(synchronize_session='fetch')

        """delete tempfile"""
        os.remove(TEMPFILE)

