"""

About

"""

__author__ = 'Anton Petrov'

import os, csv, pdb, sys, getopt, logging, smtplib, ConfigParser, collections
from email.mime.text import MIMEText
from datetime import date


from aMLSqlAlchemyClasses import *
from aLoopLoader import LoopLoader
from aPdbInfoLoader import PdbInfoLoader


class MotifAtlasLoader():
    """
    """
    def __init__(self):
        self.config = collections.defaultdict(dict)
        self.logfile = '/Users/anton/FR3D/motifatlas.log'
        self.configfile = '/Users/anton/Dropbox/Code/PyMotifLoader/motifatlas.cfg'
#         self.recalculate['distance'] = False
#         self.recalculate['coordinates'] = False
        self.import_config()

    def import_config(self):
        """
        """
        try:
            logging.info('Inside import_config')
            config = ConfigParser.RawConfigParser()
            config.read(self.configfile)
            self.config['email']['from']     = config.get('Email', 'from')
            self.config['email']['to']       = config.get('Email', 'to')
            self.config['email']['login']    = config.get('Email', 'login')
            self.config['email']['password'] = config.get('Email', 'password')
            logging.info('Leaving import_config')
            logging.info('%s', '+'*40)
        except:
            e = sys.exc_info()[1]
            logging.warning('import_config CRASHED')
            self.crash(e)


    def send_report(self):
        """
        """
        fp = open(self.logfile, 'rb')
        msg = MIMEText(fp.read())
        fp.close()

        msg['Subject'] = 'Motif Atlas Update ' + date.today().isoformat()
        server = smtplib.SMTP('smtp.gmail.com:587')
        server.ehlo()
        server.starttls()
        server.ehlo()
        server.login(self.config['email']['login'],self.config['email']['password'])
        server.sendmail(self.config['email']['from'], self.config['email']['to'], msg.as_string())
        server.quit()


    def crash(self, msg):
        logging.warning(msg)
        logging.warning('UPDATE FAILED')
        session.rollback()
        self.send_report()
        sys.exit(2)











# compute new non-redundant lists, import into the database

# extract loops, do loop QA, import into the database

# get a list of loops to be clustered

# cluster motifs, import into the database

# annotate all pdb files with these clusters

# on failure: stop, email

# log from matlab, log from python

# log with dates, clear filenames



def usage():
    print __doc__


def main(argv):
    """
    """
    logging.basicConfig(level=logging.DEBUG)
#     logging.basicConfig(filename='motifatlas.log', filemode='w', level=logging.DEBUG)
    logging.info('Initializing update')


#     try:
#         opts, args = getopt.getopt(argv, "", ['help'])
#     except getopt.GetoptError:
#         usage()
#         sys.exit(2)
#     for opt, arg in opts:
#         pass
#         if   opt == '-d':
#             M.import_distances()
# #         elif opt == '-q':
# #             L.import_loop_qa()
# #         elif opt == '-c':
# #             L.import_coordinates()
# #         elif opt == '-l':
# #             L.import_all_loops()
# #         elif opt == '-m':
# #             L.import_loop_modifications()
# #         elif opt == '-r':
# #             L.remove_loop_release(arg)
#         elif opt in ('-h', '--help'):
#             usage()
#             sys.exit()


    # mlab._dont_proxy["cell"] = True

    pdbs = ['1EKA','1HLX']#,'1S72','2AVY']


    M = MotifAtlasLoader()
    L = LoopLoader()

    """get new pdb files, import descriptions into the database"""
    P = PdbInfoLoader()
    P.update_rna_containing_pdbs()
    """import coordinates and distances into the database"""
    L.matlab_import_distances(pdbs,False)
    L.matlab_import_coordinates(pdbs,False)


    logging.info('SUCCESSFUL UPDATE')
    M.send_report()

if __name__ == "__main__":
    main(sys.argv[1:])
