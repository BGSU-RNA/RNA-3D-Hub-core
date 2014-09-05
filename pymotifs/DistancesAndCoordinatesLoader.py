"""

Loads distances and coordinates into the database.

"""

__author__ = 'Anton Petrov'

import os, csv, pdb, sys, getopt, logging
from models import session, Distances, Coordinates
from MotifAtlasBaseClass import MotifAtlasBaseClass

logger = logging.getLogger(__name__)


class DistancesAndCoordinatesLoader(MotifAtlasBaseClass):
    """
    """
    def __init__(self):
        MotifAtlasBaseClass.__init__(self)

    def import_distances(self, pdbs, recalculate=None):
        """Determines what files need to be analyzed, deletes stored data if
           necessary, loops over the pdbs, runs matlab on each of them
           independently, matlab generates a temporary csv file, it's imported
           and immediately deleted."""
        try:
            logger.info('Inside import_distances')
            if recalculate is None:
                recalculate = self.config['recalculate']['distances']
            if recalculate:
                pdb_list = pdbs
                self.__delete_distances(pdbs)
            else:
                pdb_list = self.filter_out_analyzed_pdbs(pdbs,'distances')

            if pdb_list:
                self._setup_matlab()

            for pdb_file in pdb_list:
                logger.info('Running matlab on %s', pdb_file)
                ifn, status, err_msg = self.mlab.loadDistances(pdb_file, nout=3)
                status = status[0][0]
                if status == 0:
                    self.__import_distances_from_csv(ifn)
                elif status == 2: # no nucleotides in the pdb file
                    logger.info('Pdb file %s has no nucleotides', pdb_file)
                else:
                    logger.warning('Matlab error code %i when analyzing %s',
                                     status, pdb_file)
                    self._crash(err_msg)

                self.mark_pdb_as_analyzed(pdb_file,'distances')

            logger.info('%s', '='*40)
        except:
            e = sys.exc_info()[1]
            self._crash(e)

    def __import_distances_from_csv(self, ifn):
        """Reads the csv file, imports all distances, deletes the file when done
           to avoid stale data and free up disk space"""
        logger.info('Importing distances')
        commit_every = 100000
        reader = csv.reader(open(ifn, 'rb'), delimiter=',', quotechar='"')
        for i,row in enumerate(reader):
            D = Distances(id1=row[0], id2=row[1], distance=row[2])
            try:
                session.add(D)
            except:
                logger.warning('Distance value updated')
                session.merge(D)
            """Since the files can be huge, it's unfeasible to store all
            objects in memory, have to commit regularly"""
            if i % commit_every == 0:
                session.commit()
        session.commit()
        os.remove(ifn)
        logger.info('Csv file successfully imported')

    def __delete_distances(self, pdb_list):
        """recalculate=True, so delete what's already in the database"""
        logger.info('Deleting existing records %s', ','.join(pdb_list))
        for pdb_file in pdb_list:
            session.query(Distances). \
                    filter(Distances.id1.like(pdb_file+'%')). \
                    delete(synchronize_session=False)

    def import_coordinates(self, pdbs, recalculate=None):
        """
        """
        try:
            logger.info('Inside import_coordinates')
            if recalculate is None:
                recalculate = self.config['recalculate']['coordinates']
            if recalculate is False:
                pdb_list = self.filter_out_analyzed_pdbs(pdbs,'coordinates')
            else:
                pdb_list = pdbs
                self.__delete_coordinates(pdbs)

            if pdb_list:
                self._setup_matlab()

            for pdb_file in pdb_list:
                logger.info('Running matlab on %s', pdb_file)
                ifn, status, err_msg = self.mlab.loadCoordinates(pdb_file,
                                                                 nout=3)
                status = status[0][0]
                if status == 0:
                    self.__import_coordinates_from_csv(ifn)
                elif status == 2: # no nucleotides in the pdb file
                    logger.info('Pdb file %s has no nucleotides', pdb_file)
                else:
                    logger.warning('Matlab error code %i when analyzing %s',
                                     status, pdb_file)
                    self._crash(err_msg)

                self.mark_pdb_as_analyzed(pdb_file,'coordinates')

            logger.info('%s', '='*40)
        except:
            e = sys.exc_info()[1]
            self._crash(e)

    def __import_coordinates_from_csv(self, ifn):
        """
        """
        logger.info('Importing coordinates')
        reader = csv.reader(open(ifn, 'rb'), delimiter=',', quotechar='"')
        for row in reader:
            C = Coordinates(id          = row[0],
                            pdb         = row[1],
                            pdb_type    = row[2],
                            model       = row[3],
                            chain       = row[4],
                            number      = row[5],
                            unit        = row[6],
                            ins_code    = row[7],
                            index       = row[8],
                            coordinates = row[9])
            try:
                session.add(C)
            except:
                logger.warning('Merging for %s', C.id)
                session.merge(C)
        session.commit()
        os.remove(ifn)
        logger.info('Csv file successfully imported')

    def __delete_coordinates(self, pdbs):
        """recalculate everything, delete what's already in the database"""
        logger.info('Deleting existing records before recalculation %s',
                     ','.join(pdbs))
        session.query(Coordinates).filter(Coordinates.pdb.in_(pdbs)). \
                                   delete(synchronize_session=False)


def usage():
    print __doc__


def main(argv):
    """
    """
    logging.basicConfig(level=logging.DEBUG)

    pdbs = ['1EKA','1HLX']

    D = DistancesAndCoordinatesLoader()
    D.import_distances(pdbs, recalculate=False)
    D.import_coordinates(pdbs, recalculate=False)


if __name__ == "__main__":
    main(sys.argv[1:])
