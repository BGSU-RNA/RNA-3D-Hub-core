"""

about

"""

__author__ = 'Anton Petrov'

import os, csv, pdb, sys, getopt, logging, smtplib, ConfigParser, collections
import gc
from email.mime.text import MIMEText
from datetime import date
from mlabwrap import mlab

from aMLSqlAlchemyClasses import *



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
            logging.info('+++++++++++++++++++++++++++++++++++')
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


    def import_distances(self, pdbs, distance_recalculate):
        """
        """
        try:
            logging.info('Inside import_distances')

            pdb_list = pdbs[:] # create copy, not a reference
            if distance_recalculate is False:
                """find what pbds have already been analyzed
                   and remove them from the list"""
                done = []
                for pdb_file in pdb_list:
                    if session.query(Distances). \
                               filter(Distances.id1.like(pdb_file+'%')). \
                               first() is not None:
                        done.append(pdb_file)
                [pdb_list.remove(x) for x in done]
                logging.info('Already in the database: ' + ','.join(done))
                logging.info('New files to analyze: ' + ','.join(pdb_list))
            else:
                """recalculate everything, delete what's already in the database"""
                logging.info('Deleting existing records before recalculation %s',
                             ','.join(pdb_list))
                for pdb_file in pdb_list:
                    session.query(Distances). \
                            filter(Distances.id1.like(pdb_file+'%')). \
                            delete(synchronize_session=False)
                session.commit()

            for pdb_file in pdb_list:
                ifn, status, err_msg = mlab.aGetDistances(pdb_file,nout=3)
                status = status[0][0]

                if status == 0:
                    logging.info('Pdb file %s analyzed successfully', pdb_file)
                    D = []
                    reader = csv.reader(open(ifn, 'rb'), delimiter=',', quotechar='"')
                    gc.disable()
                    for i,row in enumerate(reader):
                        session.add(Distances(id1      = row[0],
                                           id2      = row[1],
                                           distance = row[2]))
                        if i % 1000 == 0:
                            session.commit()

#                     session.add_all(D)
                    session.commit()
                    os.remove(ifn)
                    gc.enable()

                elif status == 2: # no nucleotides in the pdb file
                    logging.info('Pdb file %s has no nucleotides', pdb_file)
                else:
                    logging.warning('Matlab error code %i when analyzing %s',pdb_file,status)
                    logging.warning('Matlab says: %s', err_msg)

            logging.info('Leaving import_distances')
            logging.info('+++++++++++++++++++++++++++++++++++')
        except:
            e = sys.exc_info()[1]
            logging.warning('import_distances CRASHED')
            self.crash(e)


    def import_coordinates(self, pdbs, coordinates_recalculate):
        """
        """
        try:
            logging.info('Inside import_coordinates')

            pdb_list = pdbs[:] # create copy, not a reference
            if coordinates_recalculate is False:
                """find what pbds have already been analyzed
                   and remove them from the list"""
                done = []
                for pdb_file in pdb_list:
                    if session.query(Coordinates). \
                               filter(Coordinates.pdb==pdb_file). \
                               first() is not None:
                        done.append(pdb_file)
                [pdb_list.remove(x) for x in done]
                logging.info('Already in the database: %s', ','.join(done))
                logging.info('New files to analyze: %s', ','.join(pdb_list))
            else:
                """recalculate everything, delete what's already in the database"""
                logging.info('Deleting existing records before recalculation ' + ','.join(pdb_list))
                for pdb_file in pdb_list:
                    session.query(Coordinates). \
                            filter(Coordinates.pdb==pdb_file). \
                            delete(synchronize_session=False)
                session.commit()

            for pdb_file in pdb_list:
                ifn, status, err_msg = mlab.aGetCoordinates(pdb_file,nout=3)
                status = status[0][0]

                if status == 0:
                    logging.info('Pdb file %s analyzed successfully', pdb_file)
                    C = []
                    reader = csv.reader(open(ifn, 'rb'), delimiter=',', quotechar='"')
                    for row in reader:
                        C.append(Coordinates(id          = row[0],
                                             pdb         = row[1],
                                             pdb_type    = row[2],
                                             model       = row[3],
                                             chain       = row[4],
                                             number      = row[5],
                                             unit        = row[6],
                                             ins_code    = row[7],
                                             coordinates = row[8]))
                    session.add_all(C)
                    session.commit()
                    os.remove(ifn)

                elif status == 2: # no nucleotides in the pdb file
                    logging.info('Pdb file %s has no nucleotides', pdb_file)
                else:
                    logging.warning('Matlab error code %i when analyzing %s',pdb_file,status)
                    logging.warning('Matlab says: %s', err_msg)

            logging.info('Leaving import_coordinates')
            logging.info('+++++++++++++++++++++++++++++++++++')
        except:
            e = sys.exc_info()[1]
            logging.warning('import_coordinates CRASHED')
            self.crash(e)





# get new pdb files, import descriptions into the database

# import coordinates and distances into the database

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
    logging.basicConfig(filename='motifatlas.log', filemode='w', level=logging.DEBUG)
    logging.info('Initializing update')


    try:
        opts, args = getopt.getopt(argv, "", ['help'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    M = MotifAtlasLoader()

    for opt, arg in opts:
        pass
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

    mlab.setup()
    # mlab._dont_proxy["cell"] = True

    pdbs = ['1S72']#,'1J5E','1S72']

    M.import_distances(pdbs,False)
    M.import_coordinates(pdbs,False)
    logging.info('SUCCESSFUL UPDATE')
#     M.send_report()

if __name__ == "__main__":
    main(sys.argv[1:])
