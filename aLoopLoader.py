"""
Loop and coordinates loader

Usage: python aMLLoopLoader [options]

Options:
  -l                    import all loops
  -m                    import loop modifications
  -q                    import loop qa
  -c                    import coordinates
  -d                    import distances
  -r release_id         remove release id from loop_qa and loop_releases tables
  -h, --help            show this help

Examples:
python aMLLoopLoader.py
"""

__author__ = 'Anton Petrov'

import os, csv, pdb, sys, getopt, logging
from aMLSqlAlchemyClasses import *


class LoopLoader:
    """
    """
    def __init__(self):
        self.matlab = False
        self.f = dict()
        self.f['all_loops']   = '/FR3D/MotifAtlas/Tables/All_loops.csv'
        self.f['loop_mod']    = '/FR3D/MotifAtlas/Tables/Loops_modifications.csv'
        self.f['loop_qa']     = '/FR3D/MotifAtlas/Tables/Loops_QA.csv'
        self.f['coordinates'] = '/FR3D/MotifAtlas/Tables/Coordinates.csv'
        self.f['distances']   = '/FR3D/MotifAtlas/Tables/Distances.csv'

    def _setup_matlab(self):
        if self.matlab:
            logging.info('Matlab already running')
            return
        logging.info('Starting up matlab')
        from mlabwrap import mlab
        mlab.setup()
        self.matlab = True
        logging.info('Matlab started')
        return mlab

    def __crash(self, msg=None):
        if msg:
            logging.warning(msg)
        logging.critical('PROGRAM %s CRASHED', __name__)
        session.rollback()
        sys.exit(2)

    """ matlab_import_distances and its supporting private functions"""
    def matlab_import_distances(self, pdbs, distance_recalculate):
        """Determines what files need to be analyzed, deletes stored data if
           necessary, loops over the pdbs, runs matlab on each of them
           independently, imports csv files one by one."""
        try:
            logging.info('Inside matlab_import_distances')
            pdb_list = pdbs[:] # create a copy, not a reference
            if distance_recalculate:
                self.__delete_distances(pdb_list)
            else:
                self.__filter_done_distances(pdb_list)

            if pdb_list:
                mlab = self._setup_matlab()

            for pdb_file in pdb_list:
                logging.info('Running matlab on %s', pdb_file)
                ifn, status, err_msg = mlab.aGetDistances(pdb_file,nout=3)
                status = status[0][0]
                if status == 0:
                    self.__import_distances_from_csv(ifn)
                elif status == 2: # no nucleotides in the pdb file
                    logging.info('Pdb file %s has no nucleotides', pdb_file)
                else:
                    logging.warning('Matlab error code %i when analyzing %s',pdb_file,status)
                    logging.warning('Matlab says: %s', err_msg)

            logging.info('Leaving matlab_import_distances')
            logging.info('%s', '+'*40)
        except:
            e = sys.exc_info()[1]
            logging.critical('matlab_import_distances CRASHED')
            self.__crash(e)

    def __import_distances_from_csv(self, ifn):
        """Reads the csv file, imports all distances, deletes the file when done
           to avoid stale data and free up disk space"""
        logging.info('Importing distances')
        commit_every = 1000
        reader = csv.reader(open(ifn, 'rb'), delimiter=',', quotechar='"')
        for i,row in enumerate(reader):
            D = Distances(id1=row[0], id2=row[1], distance=row[2])
            try:
                session.add(D)
            except:
                logging.warning('Distance value updated')
                session.merge(D)
            """Since the files can be huge, it's unfeasible to store all
            objects in memory, have to commit regularly"""
            if i % commit_every == 0:
                session.commit()
        session.commit()
        os.remove(ifn)
        logging.info('Csv file successfully imported')

    def __delete_distances(self, pdb_list):
        """recalculate everything, delete what's already in the database"""
        logging.info('Deleting existing records before recalculation %s',
                     ','.join(pdb_list))
        for pdb_file in pdb_list:
            session.query(Distances). \
                    filter(Distances.id1.like(pdb_file+'%')). \
                    delete(synchronize_session=False)
        session.commit()

    def __filter_done_distances(self, pdb_list):
        """find what pbds have already been analyzed, remove them from the list"""
        logging.info('Filtering out processed distances')
        done = []
        for pdb_file in pdb_list:
            if session.query(Distances). \
                       filter(Distances.id1.like(pdb_file+'%')). \
                       first() is not None:
                done.append(pdb_file)
        [pdb_list.remove(x) for x in done]
        if done:
            logging.info('Already in the database: ' + ','.join(done))
        else:
            logging.info('No stored distances found')
        if pdb_list:
            logging.info('New files to analyze: ' + ','.join(pdb_list))
        else:
            logging.info('No new files to analyze')

    """matlab_import_coordinates and its supporting private functions"""
    def matlab_import_coordinates(self, pdbs, coordinates_recalculate):
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
            logging.warning('import_coordinates __crashED')
            self.__crash(e)

    def import_coordinates(self):
        ifn = self.f['coordinates']
        reader = csv.reader(open(ifn, 'rb'), delimiter=',', quotechar='"')
        for row in reader:
            coord = Coordinates(parts = row)
            c = session.query(Coordinates).filter(Coordinates.nt_id == coord.nt_id).first()
            if not c:
                session.add(coord)
            elif c.coordinates != coord.coordinates:
                c.coordinates = coord.coordinates
#             else:
#                 print 'Residue ', coord.nt_id, ' already imported'

        session.commit()
        print 'Coordinates imported'

    def import_distances(self):
        ifn = self.f['distances']
        reader = csv.reader(open(ifn, 'rb'), delimiter=',', quotechar='"')
        for row in reader:
            dist = Distances(parts = row)
            d = session.query(Distances).filter(Distances.id1 == d.id1) \
                                        .filter(Distances.id2 == d.id2).first()
            if not d:
                session.add(dist)
            elif d.distance != dist.distance:
                d.distance = dist.distance
#             else:
#                 print 'Distance between %s and %s already imported', %(dist.id1, dist.id2)

        session.commit()
        print 'Distances imported'


    def import_all_loops(self):
        """
        """
        ifn = self.f['all_loops']
        reader = csv.reader(open(ifn, 'rb'), delimiter=',', quotechar='"')
        for row in reader:
            loop = AllLoops(parts = row)
            l = session.query(AllLoops).filter(AllLoops.id == loop.id).first()
            if not l:
                session.add(loop)
            elif l != loop:
                updated += 1
                l = loop

        session.commit()
        print 'Loops imported'


    def import_loop_modifications(self):
        """
        """
        ifn = self.f['loop_mod']
        reader  = csv.reader(open(ifn, 'rb'), delimiter=',', quotechar='"')
        updated = unchanged = new = 0
        try:
            for row in reader:
                mod = LoopModifications(id=row[0], modification=row[1])
                m = session.query(LoopModifications).filter(LoopModifications.id == mod.id).first()

                if not m:
                    new += 1
                    session.add(mod)
                elif m.modification != mod.modification:
                    print 'Old modification: ', m.modification
                    print 'New modification: ', mod.modification
                    m.modification = mod.modification
                    updated += 1
                else:
                    unchanged += 1
        except csv.Error, e:
            print 'file %s, line %d: %s' % (ifn, reader.line_num, e)

        session.commit()
        print 'Loop modifications imported, %i updated, %i unchanged, %i new' % (updated, unchanged, new)


    def import_loop_qa(self):
        """
        """
        ifn = self.f['loop_qa']
        reader = csv.reader(open(ifn, 'rb'))
        imported = updated = 0
        release = LoopRelease(mode='minor',description='test')
        try:
            for row in reader:
                qa = LoopQA(id=row[0], code=row[1], release_id=release.id)
                q = session.query(LoopQA).filter(LoopQA.id == qa.id) \
                                         .filter(LoopQA.release_id == release.id).first()
                if not q:
                    session.add(qa)
                    imported += 1
                else:
                    updated += 1
                    if row[1] == '1':
                        q.valid       = 1
                    elif row[1] == '2':
                        q.missing_nt  = 1
                    elif row[1] == '3':
                        q.modified_nt = 1
                    elif row[1] == '4':
                        q.complementary = 1


        except csv.Error, e:
            print 'file %s, line %d: %s' % (ifn, reader.line_num, e)

        session.add(release)
        session.commit()
        print 'Loop QA imported, %i loops in release %s, %i updated' % (imported, release.id, updated)

    def remove_loop_release(self, id):
        """
        """
        print 'Are you sure you want to delete release %s?' % id
        print 'Press "c" to confirm, "q" to quit'
        pdb.set_trace()
        session.query(LoopQA).filter(LoopQA.release_id == id).delete()
        session.query(LoopRelease).filter(LoopRelease.id == id).delete()
        logging.info('Release %s successfully removed', id)




def usage():
    print __doc__

def main(argv):
    """
    """
    try:
        opts, args = getopt.getopt(argv, "qlcdmhr:", ['help'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

#     logging.basicConfig(filename='motifatlas.log', filemode='w', level=logging.DEBUG)
    logging.basicConfig(level=logging.DEBUG)
    L = LoopLoader()

#     mlab.setup()
    pdbs = ['1S72','124D','2AW4']

    for opt, arg in opts:
        if   opt == '-d':
#             L.import_distances()
            L.matlab_import_distances(pdbs,False)
        elif opt == '-q':
            L.import_loop_qa()
        elif opt == '-c':
            L.import_coordinates()
        elif opt == '-l':
            L.import_all_loops()
        elif opt == '-m':
            L.import_loop_modifications()
        elif opt == '-r':
            L.remove_loop_release(arg)
        elif opt in ('-h', '--help'):
            usage()
            sys.exit()
        else:
            usage()
            sys.exit()


if __name__ == "__main__":
    main(sys.argv[1:])