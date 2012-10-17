"""



"""

__author__ = 'Anton Petrov'

import sys
import logging
import smtplib
import ConfigParser
import collections
import datetime
import os
from email.mime.text import MIMEText
from time import localtime, strftime

from models import session, PdbAnalysisStatus


class MotifAtlasBaseClass:
    """
        Don't use logging anywhere in this constructor or the functions it calls
    """
    def __init__(self):
        self.mlab   = False
        self.config = collections.defaultdict(dict)
        script_path = os.path.dirname(os.path.abspath( __file__ ))
        self.configfile = os.path.join(script_path, 'motifatlas.cfg')
        self.import_config()
        self.log = ''
        self.log_filename = 'rna3dhub_log.txt'

    def start_logging(self):
        """
            Overwrites the old log file.
        """
        log_dir = self.config['locations']['log_dir']
        self.log = os.path.join(log_dir, self.log_filename)
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)
        logging.basicConfig(filename=self.log,
                            level=logging.DEBUG,
                            filemode='w')
        print 'Log file %s' % self.log

    def _setup_matlab(self):
        if self.mlab:
            return
        logging.info('Starting up matlab')
        os.chdir(self.config['locations']['fr3d_root'])
        from mlabwrap import mlab
        self.mlab = mlab
        self.mlab.setup() # add matlab paths
        # self.mlab._dont_proxy["cell"] = True
        logging.info('Matlab started')

    def filter_out_analyzed_pdbs(self, pdbs, column_name):
        """Checks whether the pdb files were processed . Returns only the files
        that need to be analyzed. The `column_name` parameter corresponds to
        the column name of the pdb_analysis_status table."""
        pdb_list = pdbs[:] # copy, not reference
        done = session.query(PdbAnalysisStatus). \
                       filter(PdbAnalysisStatus.id.in_(pdbs)). \
                       filter(getattr(PdbAnalysisStatus,column_name.lower()) != None). \
                       all()
        [pdb_list.remove(x.id) for x in done]
        if pdb_list:
            logging.info('New files to analyze: ' + ','.join(pdb_list))
        else:
            logging.info('No new files to analyze')
        return pdb_list

    def mark_pdb_as_analyzed(self, pdb_id, column_name):
        """
        """
        P = PdbAnalysisStatus(id = pdb_id)
        setattr(P,column_name.lower(),datetime.datetime.now())
        session.merge(P)
        session.commit()
        logging.info('Updated %s status for pdb %s', column_name, pdb_id)

    def import_config(self):
        """
        """
        try:
            config = ConfigParser.RawConfigParser()
            config.read(self.configfile)
            """general settings"""
            section = 'general'
            keys = ['environment']
            for k in keys: self.config[section][k] = config.get(section,k)
            """email settings"""
            section = 'email'
            keys = ['from','to','login','password','subject']
            for k in keys: self.config[section][k] = config.get(section,k)
            """recalculation settings"""
            section = 'recalculate'
            keys = ['coordinates','distances','interactions','IL','HL','J3',
                    'redundant_nts','best_chains_and_models']
            for k in keys: self.config[section][k] = config.getboolean(section,k)
            """logging"""
            self.config['logfile'] = 'motifatlas.log'
            """locations"""
            section = 'locations'
            keys = ['loops_mat_files', 'loops_search_dir', 'log_dir',
                    'releases_dir', 'nrlists_dir', 'fr3d_root',
                    '2ds_destination', 'mlab_app']
            for k in keys: self.config[section][k] = config.get(section,k)
            """release modes"""
            section = 'release_mode'
            keys = ['loops','motifs','nrlist']
            for k in keys: self.config[section][k] = config.get(section,k)
        except:
            e = sys.exc_info()[1]
            self._crash(e)

    def send_report(self):
        """
        """
        try:
            fp = open(self.log, 'rb')
            msg = MIMEText(fp.read())
            fp.close()
            msg['Subject'] = ' '.join([self.config['email']['subject'],
                                 strftime("%Y-%m-%d", localtime())])
            server = smtplib.SMTP('smtp.gmail.com:587')
            server.ehlo()
            server.starttls()
            server.ehlo()
            server.login(self.config['email']['login'],
                         self.config['email']['password'])
            server.sendmail(self.config['email']['from'],
                            self.config['email']['to'],
                            msg.as_string())
            server.quit()
        except:
            e = sys.exc_info()[1]
            logging.critical(e)
            sys.exit(2)

    def _crash(self, msg=None):
        """
        """
        if msg:
            logging.critical(msg)
        try:
            session.rollback()
        except:
            pass
        self.send_report()
        sys.exit(2)
