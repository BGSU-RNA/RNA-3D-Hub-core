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

import os, csv, pdb, sys, getopt
from aMLSqlAlchemyClasses import *


class LoopLoader:
    """
    """
    def __init__(self):
        self.f = dict()
        self.f['all_loops']   = '/FR3D/MotifAtlas/Tables/All_loops.csv'
        self.f['loop_mod']    = '/FR3D/MotifAtlas/Tables/Loops_modifications.csv'
        self.f['loop_qa']     = '/FR3D/MotifAtlas/Tables/Loops_QA.csv'
        self.f['coordinates'] = '/FR3D/MotifAtlas/Tables/Coordinates.csv'
        self.f['distances']   = '/FR3D/MotifAtlas/Tables/Distances.csv'


    def remove_loop_release(self, id):
        """
        """
        print 'Are you sure you want to delete release %s? Press c to confirm, q to quit' % id
        pdb.set_trace()
        session.query(LoopQA).filter(LoopQA.release_id == id).delete()
        session.query(LoopRelease).filter(LoopRelease.id == id).delete()
        print 'Release %s successfully removed' % id


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

    L = LoopLoader()

    for opt, arg in opts:
        if   opt == '-d':
            L.import_distances()
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


if __name__ == "__main__":
    main(sys.argv[1:])
