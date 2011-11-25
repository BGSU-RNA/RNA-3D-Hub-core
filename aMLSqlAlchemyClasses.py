"""

"""

import random, datetime, math, sys, pdb, csv, os, shutil, collections

from sqlalchemy import *
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

import sqlalchemy.exc

engine  = create_engine('mysql://root:bioinfo@localhost/mltest')
Session = sessionmaker(bind=engine)
session = Session()

Base = declarative_base()


################################################################################
# Auxiliary functions
################################################################################

def drop_all():
    Base.metadata.drop_all(engine)

def list_done():
    done = []
    for release in session.query(Release).all():
        done.append(release.description)
    return done

def list_all_releases():
    all_releases = []
    for release in session.query(Release).order_by(desc(Release.date)).all():
        all_releases.append(release.id)
    return all_releases

################################################################################
# Motif tables declarations
################################################################################

class AllLoops(Base):
    """
    """
    __tablename__ = 'loops_all'
    # "IL_157D_001","IL","157D","001","6","CGA*UAG","UAG*CGA","G*A","A*G","157D.pdb","fullid"
    id            = Column(String(11), primary_key=True)
    type          = Column(String(2))
    pdb           = Column(String(4))
    sequential_id = Column(String(3))
    length        = Column(Integer)
    seq           = Column(Text)
    r_seq         = Column(Text)
    nwc_seq       = Column(Text)
    r_nwc_seq     = Column(Text)
    pdb_file      = Column(String(10))
    nt_ids        = Column(Text)
    loop_name     = Column(Text)

    def __init__(self, parts):
        self.id   = parts[0]
        self.type = parts[1]
        self.pdb  = parts[2]
        self.sequential_id = parts[3]
        self.length        = parts[4]
        self.seq           = parts[5]
        self.r_seq     = parts[6]
        self.nwc_seq   = parts[7]
        self.r_nwc_seq = parts[8]
        self.pdb_file  = parts[9]
        self.nt_ids    = parts[10]
        self.loop_name = parts[11]

    def __repr__(self):
        return "<Loop('%s','%s')>" % (self.id, self.seq)


class LoopModifications(Base):
    """
    """
    __tablename__ = 'loop_modifications'
    # IL_1C2W_062,7MG
    id           = Column(String(11), primary_key=True)
    modification = Column(String(3))

    def __init__(self, id = '', modification = ''):
        self.id = id
        self.modification = modification

    def __repr__(self):
        return "<Modified('%s','%s')>" % (self.id, self.modification)


class LoopQA(Base):
    """
    """
    __tablename__ = 'loop_qa'
    # IL_1C2W_062,1
    # 1 - valid, 2 - missing, 3 - modified
    id    = Column(String(11), primary_key=True)
    valid = Column(Boolean)
    modified_nt    = Column(Boolean)
    missing_nt     = Column(Boolean)
    complementary  = Column(Boolean);
    release_id     = Column(String(4), primary_key=True)

    def __init__(self, id='', code='0', release_id=''):
        self.id = id
        self.release_id = release_id
        if code == '1':
            self.valid         = 1
            self.modified_nt   = 0
            self.missing_nt    = 0
            self.complementary = 0
        elif code == '2':
            self.valid         = 0
            self.modified_nt   = 0
            self.missing_nt    = 1
            self.complementary = 0
        elif code == '3':
            self.valid         = 0
            self.modified_nt   = 1
            self.missing_nt    = 0
            self.complementary = 0
        elif code == '4':
            self.valid         = 0
            self.modified_nt   = 0
            self.missing_nt    = 0
            self.complementary = 1
        else:
            self.valid         = 0
            self.modified_nt   = 0
            self.missing_nt    = 0
            self.complementary = 0

    def __repr__(self):
        return "<LoopQA('%s','%s','%s')>" % (self.id, self.valid,self.release_id)


class LoopRelease(Base):
    """
    """
    __tablename__ = 'loop_releases'

    id          = Column(String(4), primary_key=True)
    date        = Column(DateTime)
    description = Column(Text)

    def __init__(self, mode='', description=''):
        self.description = description
        self.mode = mode
        self.compute_new_release_id()
        self.get_date()

    def __repr__(self):
        return "<LoopRelease('%s','%s','%s')>" % (self.id, self.date, self.description)

    def compute_new_release_id(self):
        prev = session.query(LoopRelease).order_by(desc(LoopRelease.date)).first()
        if prev is None:
            self.id = 0.1
        elif self.mode == 'major':
            parts = prev.id.split('.')
            self.id = '.'.join([str(int(parts[0])+1), '0'])
        else:
            parts = prev.id.split('.')
            self.id = '.'.join([parts[0], str(int(parts[1])+1)])

    def get_date(self):
        self.date = datetime.datetime.now()


class Release(Base):
    """
    """
    __tablename__ = 'ml_releases'

    id          = Column(String(4), primary_key=True)
    date        = Column(DateTime)
    description = Column(Text)

    def __init__(self, mode='', description=''):
        self.description = description
        self.mode = mode
        self.compute_new_release_id()
        self.get_date()

    def __repr__(self):
        return "<Release('%s','%s','%s')>" % (self.id, self.date, self.description)

    def compute_new_release_id(self):
        prev = session.query(Release).order_by(desc(Release.date)).first()
        if prev is None:
            self.id = 0.1
        elif self.mode == 'major':
            parts = prev.id.split('.')
            self.id = '.'.join([str(int(parts[0])+1), '0'])
        else:
            parts = prev.id.split('.')
            self.id = '.'.join([parts[0], str(int(parts[1])+1)])

    def get_date(self):
        self.date = datetime.datetime.now()


class Parents(Base):
    __tablename__ = 'ml_history'

    motif_id   = Column(String(11), primary_key=True)
    release_id = Column(String(4),  primary_key=True )
    parents    = Column(Text)

    def __repr__(self):
        return "<Parents('%s','%s','%s')>" % (self.motif_id, self.release_id, self.parents)


class Loop(Base):
    __tablename__ = 'ml_loops'

    id         = Column(String(11), primary_key=True)
    motif_id   = Column(String(11), primary_key=True)
    release_id = Column(String(4),  primary_key=True)

    def __repr__(self):
        return "<Loop('%s','%s','%s')>" % (self.id, self.motif_id, self.release_id)


class Motif(Base):
    """
    """
    __tablename__ = 'ml_motifs'

    id         = Column(String(11), primary_key=True)
    release_id = Column(String(4),  primary_key=True)
    type       = Column(String(2)) # IL, HL, JL
    handle     = Column(String(5)) # XXXXX
    version    = Column(Integer)

    def __init__(self, id='', release_id='', increment=False):
        self.release_id = release_id
        if id == '':
            self.get_new_motif_id()
        elif increment is True:
            self.increment_motif_id(id)
        else:
            self.populate_fields(id)

    def __repr__(self):
        return "<Motif('%s','%s','%s')>" % (self.id, self.type, self.release_id)

    def get_new_motif_id(self):
        self.type = 'IL'
        while True:
            self.handle = '%05d' % random.randrange(99999)
            if session.query(Motif).filter(Motif.handle==self.handle).first() is None:
                break
        self.version = 1
        self.id = self.type + '_' + self.handle + '.' + str(self.version)

    def increment_motif_id(self, id):
        self.type = 'IL'
        self.handle  = id[3:8]
        self.version = int(id[9:]) + 1
        self.id = self.type + '_' + self.handle + '.' + str(self.version)

    def populate_fields(self, id):
        self.id = id
        self.type = 'IL'
        self.handle  = id[3:8]
        self.version = int(id[9:])



class LoopAnnotation(Base):
    __tablename__ = 'ml_loop_annotations'

    id          = Column(String(11), primary_key=True)
    common_name = Column(Text)
    reference   = Column(Text)


class SetDiff(Base):
    __tablename__ = 'ml_set_diff'

    motif_id1  = Column(String(11), primary_key=True)
    motif_id2  = Column(String(11), primary_key=True)
    release_id = Column(String(4),  primary_key=True)
    intersection  = Column(Text)
    overlap       = Column(Float)
    one_minus_two = Column(Text)
    two_minus_one = Column(Text)

    def __repr__(self):
        return "<SetDiff('%s','%s','%s')>" % (self.motif_id1, self.motif_id2, self.release_id)


class LoopOrder(Base):
    __tablename__ = 'ml_loop_order'

    motif_id         = Column(String(11), primary_key=True)
    loop_id          = Column(String(11), primary_key=True)
    release_id       = Column(String(4),  primary_key=True)
    original_order   = Column(Integer)
    similarity_order = Column(Integer)

    def __repr(self):
        return "<LoopOrder('%s','%s','%i')>" % (self.motif_id, self.loop_id, self.original_order)


class LoopPosition(Base):
    __tablename__ = 'ml_loop_positions'

    motif_id   = Column(String(11), primary_key=True)
    loop_id    = Column(String(11), primary_key=True)
    release_id = Column(String(4),  primary_key=True)
    nt_id      = Column(String(30), primary_key=True)
    position   = Column(Integer)

    def __repr__(self):
        return "<LoopPosition('%s','%s','%i')>" % (self.loop_id, self.nt_id, self.position)


class LoopDiscrepancy(Base):
    __tablename__ = 'ml_mutual_discrepancy'

    loop_id1    = Column(String(11), primary_key=True)
    loop_id2    = Column(String(11), primary_key=True)
    release_id  = Column(String(4),  primary_key=True)
    discrepancy = Column(Float)

    def __repr__(self):
        return "<LoopDiscrepancy('%s','%s','%f')>" % (self.loop_id1, self.loop_id2, self.discrepancy)

################################################################################
# End of motif tables declarations
################################################################################

################################################################################
# Collection comparison classes
################################################################################

class MotifCollection:
    """

    """
    def __init__(self, file='', release=''):
        self.loops  = []
        self.groups = []
        self.sl     = set()
        self.sg     = set()
        self.d      = collections.defaultdict(list) # [group: {loopX..loopY} ]
        self.sd     = collections.defaultdict(set)  # [group: set({loopX..loopY}]

        if file is not '':
            self.file = file
            self.read_motifs_csv_file()
        elif release == 'latest':
            self.get_latest_release()
        else:
            self.release_id = release
            self.get_release()

        self.make_dictionary()
        self.make_sets()

    def read_motifs_csv_file(self):
        r = csv.reader(open(self.file), delimiter=',', quotechar='"')
        for row in r:
            self.loops.append(row[0])
            self.groups.append(row[1])

    def make_sets(self):
        self.sl = set(self.loops)
        self.sg = set(self.groups)
        for k,v in self.d.iteritems():
            self.sd[k] = set(v)

    def make_dictionary(self):
        for i,loop in enumerate(self.loops):
            self.d[self.groups[i]].append(loop)

    def get_latest_release(self):
        if session.query(Release).first() is None:
            print 'No previous releases found'
            return
        release = session.query(Release).order_by(desc(Release.date)).first()
        loops = session.query(Loop).filter(Loop.release_id==release.id)
        for loop in loops:
            self.loops.append(loop.id)
            self.groups.append(loop.motif_id)

    def get_release(self):
        for loop in session.query(Loop).filter(Loop.release_id==self.release_id).all():
            self.loops.append(loop.id)
            self.groups.append(loop.motif_id)



class MotifCollectionMerger:
    """
        Analyzes two collections and establishes correspondences.
        c1 = Collection1, new collection
        c2 = Collection2, old collection
    """

    def __init__(self, c1, c2):
        self.c1 = c1
        self.c2 = c2

        self.intersection   = collections.defaultdict(dict)
        self.noIntersection = dict()
        self.overlap        = collections.defaultdict(dict)
        self.setdiff        = collections.defaultdict(dict)

        self.correspond  = dict()
        self.exact_match = dict()
        self.new_ids     = []
        self.parents     = collections.defaultdict(list)
        self.explanation = dict()

        self.minOverlap = float(2)/3

        self.compare_releases()
        self.show_report()
        self.establish_correspondences()


    def compare_releases(self):
        """
        """
        for groupId, group in self.c1.sd.iteritems():
            for motifId, motif in self.c2.sd.iteritems():
                t = group.intersection(motif)
                if len(t) > 0:
                    self.intersection[groupId][motifId] = t
                    self.intersection[motifId][groupId] = t
                    self.setdiff[groupId][motifId] = group - motif
                    self.setdiff[motifId][groupId] = motif - group
                    self.overlap[groupId][motifId] = float(len(t)) / len(group)
                    self.overlap[motifId][groupId] = float(len(t)) / len(motif)
                else:
                    self.noIntersection[groupId] = motifId
                    self.noIntersection[motifId] = groupId


    def establish_correspondences(self):
        """
        assign new id if:
            more than 2 parents
            intersection < 2/3
            the size of the final cluster is more than 2/3 larger
        """
        for groupId in self.c1.sg:
            # if intersected...
            if groupId in self.intersection:
                match = self.intersection[groupId]
                # with exactly one motif
                if len(match) == 1:
                    motifId = match.keys()[0]
                    # and they both didn't match anything else - assign the same id
                    if self.overlap[groupId][motifId] == 1 and \
                       self.overlap[motifId][groupId] == 1:
                       self.exact_match[groupId] = motifId
                       self.explanation[groupId] = 'Exact match'
                       print 'Groups %s and motif %s match up exactly' % (groupId, motifId)
                    # and the match is not perfect
                    else:
                        self.parents[groupId].append(motifId)
                        if self.overlap[groupId][motifId] >= self.minOverlap and \
                           self.overlap[motifId][groupId] >= self.minOverlap:
                            self.correspond[groupId] = motifId
                            self.explanation[groupId] = 'Updated, 1 parent'
                        else:
                            self.new_ids.append(groupId)
                            self.parents[groupId].append(motifId)
                            self.explanation[groupId] = 'New id, 1 parent'

                # with two motifs
                elif len(match) == 2:
                    success = False
                    for motifId,v in match.iteritems():
                        self.parents[groupId].append(motifId)
                        # if the overlap is big enough - assign the same id
                        if self.overlap[groupId][motifId] >= self.minOverlap and \
                           self.overlap[motifId][groupId] >= self.minOverlap:
                            success = True
                            self.correspond[groupId] = motifId
                            self.explanation[groupId] = 'Updated, 2 parents'
                            break
                    # if the overlap is not big enough - assign new id
                    if success == False:
                        self.new_ids.append(groupId)
                        self.explanation[groupId] = 'New id, 2 parents'
                    else:
                        # don't need to do anything because success==True
                        pass

                # with more than 2 motifs - assign new id
                elif len(match) > 2:
                    print 'Group %s has more than 2 parents' % groupId
                    self.explanation[groupId] = '> 2 parents'
                    self.new_ids.append(groupId)
                    for motifId,v in match.iteritems():
                        self.parents[groupId].append(motifId)


            # if didn't intersect with anything
            else:
                self.explanation[groupId] = 'New id, no parents'
                self.new_ids.append(groupId)


    def show_report(self):
        """
        """
        print len(self.c1.loops), ' new loops'
        print len(self.c2.loops), ' old loops'
        print len(self.c1.d.keys()), ' new groups'
        print len(self.c2.d.keys()), ' old groups'

        print len(self.c1.sl - self.c2.sl), ' new loops not in the database'
        print len(self.c2.sl - self.c1.sl), ' old loops not in the new list'

        all_new = set(self.c1.groups) - set(self.intersection.keys())
        if len(all_new) > 0:
            print len(all_new), ' entirely new groups'
            print [x for x in all_new]

        only_old = set(self.c2.groups) - set(self.intersection.keys())
        if len(only_old) > 0:
            print len(only_old), ' motifs only in the database'
            print [x for x in only_old]

################################################################################
# End of collections comparison classes
################################################################################


################################################################################
# Database operations
################################################################################

class Uploader:
    """
    """
    def __init__(self, ensembles=None, mode='', description='',files={}):
        self.c         = ensembles
        self.release   = Release(mode=mode, description=description)
        self.motifs    = []
        self.loops     = []
        self.history   = []
        self.intersection = []
        self.final_ids = dict()
        self.files = files

        self.loop_order       = []
        self.loop_positions   = []
        self.loop_discrepancy = []

        self.__finalize()
        self.__process_set_diff()
        self.__process_motif_loop_order()
        self.__process_motif_loop_positions()
        self.__process_mutual_discrepancy()
        self.show_report()
        self.__commit()
        self.__rename_mat_files()


    def __finalize(self):
        """
        """
        for group_id in self.c.c1.sg:
            if group_id in self.c.new_ids:
                motif = Motif(release_id=self.release.id)
                if self.c.parents.has_key(group_id):
                    parents = ','.join(self.c.parents[group_id])
                else:
                    parents = ''
                print 'Group %s assigned new id %s' % (group_id, motif.id)
            elif group_id in self.c.correspond:
                old_id  = self.c.correspond[group_id]
                motif   = Motif(id=old_id, release_id=self.release.id, increment=True)
                parents = ','.join(set([old_id] + self.c.parents[group_id]))
                print 'Group %s corresponds to motif %s and is assigned new id %s' % (group_id, old_id, motif.id)
            elif group_id in self.c.exact_match:
                id = self.c.exact_match[group_id]
                motif = Motif(id=id, release_id=self.release.id)
                parents = ''
                print 'Group %s matches exactly motif %s' % (group_id, motif.id)
            else:
                print 'Major problem'

            self.motifs.append(motif)
            if parents != '':
                self.history.append(Parents(motif_id=motif.id, release_id=self.release.id, parents=parents))
            self.final_ids[group_id] = motif.id
            for loop_id in self.c.c1.d[group_id]:
                self.loops.append(Loop(id=loop_id, motif_id=motif.id, release_id=self.release.id))


    def __process_set_diff(self):
        """
        """
        for loop_id in self.c.c1.sg:
            for motif_id in self.c.c2.sg:
                if motif_id != self.final_ids[loop_id] and \
                   self.c.intersection.has_key(loop_id) and \
                   self.c.intersection[loop_id].has_key(motif_id):

                    self.intersection.append( SetDiff(
                                motif_id1 = self.final_ids[loop_id],
                                motif_id2 = motif_id,
                                release_id = self.release.id,
                                intersection = ','.join(self.c.intersection[loop_id][motif_id]),
                                overlap = self.c.overlap[loop_id][motif_id],
                                one_minus_two = ','.join(self.c.setdiff[loop_id][motif_id]),
                                two_minus_one = ','.join(self.c.setdiff[motif_id][loop_id])
                    ))
                    self.intersection.append( SetDiff(
                                motif_id1 = motif_id,
                                motif_id2 = self.final_ids[loop_id],
                                release_id = self.release.id,
                                intersection = ','.join(self.c.intersection[loop_id][motif_id]),
                                overlap = self.c.overlap[motif_id][loop_id],
                                one_minus_two = ','.join(self.c.setdiff[motif_id][loop_id]),
                                two_minus_one = ','.join(self.c.setdiff[loop_id][motif_id])
                    ))



    def show_report(self):
        print 'Release\n',      self.release,      '\n======================='
        print 'Motifs\n',       self.motifs,       '\n======================='
        print 'Loops\n',        self.loops,        '\n======================='
        print 'History\n',      self.history,      '\n======================='
        print 'Intersection\n', self.intersection, '\n======================='

    def remove_release(self, release):
        session.query(Release).filter(Release.id==release).delete()
        session.query(Motif).filter(Motif.release_id==release).delete()
        session.query(Loop).filter(Loop.release_id==release).delete()
        session.query(SetDiff).filter(SetDiff.release_id==release).delete()
        session.query(Parents).filter(Parents.release_id==release).delete()
        session.query(LoopDiscrepancy).filter(LoopDiscrepancy.release_id==release).delete()

    def __commit(self):
        failed = False
        try:
            session.add(self.release)
            session.add_all(self.motifs)
            session.add_all(self.loops)
            session.add_all(self.history)
            session.add_all(self.intersection)

            session.add_all(self.loop_order)
            session.add_all(self.loop_positions)
            session.add_all(self.loop_discrepancy)
            session.commit()
            print 'Successful update'
        except sqlalchemy.exc.SQLAlchemyError, e:
            print 'Update failed. Rolling back.'
            print str(e)
            # print sys.exc_info()[0]
            session.rollback()
            self.remove_release(self.release.id)
        except sqlalchemy.exc.DBAPIError, e:
            print 'Update failed. Rolling back.'
            print str(e)
            # print sys.exc_info()[0]
            session.rollback()
            self.remove_release(self.release.id)
        except sys.exc_info()[0]:
            print 'Update failed. Rolling back.'
            print sys.exc_info()[0]
            session.rollback()
            self.remove_release(self.release.id)

    def __process_motif_loop_order(self):
        """
        IL_017,IL_1O9M_003,1,1
        """
        r = csv.reader(open(self.files['MotifLoopOrder']), delimiter=',', quotechar='"')
        for row in r:
            self.loop_order.append(LoopOrder(motif_id=self.final_ids[row[0]],
                                        loop_id=row[1],
                                        release_id=self.release.id,
                                        original_order=row[2],
                                        similarity_order=row[3]
                                   ))


    def __process_motif_loop_positions(self):
        """
        IL_017,IL_1O9M_003,1O9M_AU_1_A_1491_G_,1
        """
        r = csv.reader(open(self.files['MotifPositions']), delimiter=',', quotechar='"')
        for row in r:
            self.loop_positions.append(LoopPosition(motif_id=self.final_ids[row[0]],
                                        loop_id=row[1],
                                        release_id=self.release.id,
                                        nt_id=row[2],
                                        position=row[3]
                                   ))


    def __process_mutual_discrepancy(self):
        """
        IL_1O9M_003,0.0000,IL_1O9M_003
        """
        r = csv.reader(open(self.files['MutualDiscrepancy']), delimiter=',', quotechar='"')
        for row in r:
            self.loop_discrepancy.append(LoopDiscrepancy(loop_id1=row[0],
                                        discrepancy=row[1],
                                        loop_id2=row[2],
                                        release_id=self.release.id,
                                   ))

    def __rename_mat_files(self):

        if not os.path.exists(self.files['MatFiles']):
            os.mkdir(self.files['MatFiles'])

        for file in self.c.c1.sg:
            src = os.path.join(self.files['folder'],file+'.mat')
            dst = os.path.join(self.files['MatFiles'], self.final_ids[file] + '.mat')
            if os.path.exists(src):
                shutil.copyfile(src, dst)
            else:
                print "File %s wasn't found" % src

################################################################################
# End of database operations
################################################################################


Base.metadata.create_all(engine)
