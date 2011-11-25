"""

"""

import random, datetime, math, sys, pdb, csv, collections

from sqlalchemy import *
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

import sqlalchemy.exc

engine  = create_engine('mysql://root:bioinfo@localhost/nrtest')
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
    for release in session.query(NR_release).all():
        done.append(release.description)
    return done

def list_all_releases():
    all_releases = []
    for release in session.query(NR_release).order_by(desc(NR_release.date)).all():
        all_releases.append(release.id)
    return all_releases


################################################################################
# NR tables declarations
################################################################################

class NR_pdb(Base): # = loop
    """
    """
    __tablename__ = 'nr_pdbs'
    id         = Column(String(4),  primary_key=True)
    class_id   = Column(String(17), primary_key=True)
    release_id = Column(String(4),  primary_key=True)
    rep        = Column(Boolean)

    def __repr__(self):
        return "<PDB('%s','%s','%s')>" % (self.id, self.class_id, self.release_id)


class NR_class(Base): # = motif
    """
    """
    __tablename__ = 'nr_classes'

    id         = Column(String(17), primary_key=True)
    release_id = Column(String(4),  primary_key=True)
    resolution = Column(String(4))
    handle     = Column(String(5))
    version    = Column(String(3))
    comment    = Column(Text)

    def __init__(self, id='', increment=False, release_id='', resolution='', comment=''):
        self.release_id = release_id
        self.resolution = resolution
        self.comment    = comment
        if id == '':
            self.get_new_motif_id()
        elif increment is True:
            self.increment_motif_id(id)
        else:
            self.populate_fields(id)

    def generate_id(self):
        self.id = '.'.join(['_'.join([self.type,self.resolution,self.handle]),str(self.version)])

    def get_new_motif_id(self):
        self.type = 'NR'
        while True:
            self.handle = '%05d' % random.randrange(99999)
            if session.query(NR_handle).filter(NR_handle.id==self.handle).first() is None:
                if session.query(NR_class).filter(NR_class.handle==self.handle).first() is None:
                    h = NR_handle(id=self.handle)
                    session.add(h)
                    break
        self.version = 1
        self.generate_id()

    def populate_fields(self, id):
        self.id = id
        [self.type, self.resolution, handle_version] = id.split('_')
        [self.handle, version] = handle_version.split('.')
        self.version = int(version)

    def increment_motif_id(self, id):
        self.populate_fields(id)
        self.version = int(self.version) + 1
        self.generate_id()

    def __repr__(self):
        return "<NR_class('%s','%s','%s')>" % (self.id, self.type, self.release_id)


class NR_handle(Base):
    """
    """
    __tablename__ = 'nr_handles'
    id = Column(String(5), primary_key=True)


class NR_release(Base):
    """
    mode = minor/major/reuse. Release mode
    """
    __tablename__ = 'nr_releases'

    id          = Column(String(4), primary_key=True)
    date        = Column(DateTime)
    description = Column(Text)

    def __init__(self, mode='minor', description=''):
        self.description = description
        self.mode = mode
        self.compute_new_release_id()
        self.get_date()

    def __repr__(self):
        return "<NR_release('%s','%s','%s')>" % (self.id, self.date, self.description)

    def compute_new_release_id(self):
        prev = session.query(NR_release).order_by(desc(NR_release.date)).first()
        if prev is None:
            self.id = '0.1'
        elif self.mode == 'major':
            parts = prev.id.split('.')
            self.id = '.'.join([str(int(parts[0])+1), '0'])
        elif self.mode == 'reuse':
            self.id = prev.id
        else:
            parts = prev.id.split('.')
            self.id = '.'.join([parts[0], str(int(parts[1])+1)])

    def get_date(self):
        self.date = datetime.datetime.now()


class NR_setdiff(Base):
    """
    """
    __tablename__ = 'nr_set_diff'

    nr_class1     = Column(String(17), primary_key=True)
    nr_class2     = Column(String(17), primary_key=True)
    release_id    = Column(String(4),  primary_key=True)
    intersection  = Column(Text)
    overlap       = Column(Float)
    one_minus_two = Column(Text)
    two_minus_one = Column(Text)

    def __repr__(self):
        return "<NRSetDiff('%s','%s','%s')>" % (self.nr_class1, self.nr_class2, self.release_id)


class NR_parents(Base):
    """
    """
    __tablename__ = 'nr_parents'

    class_id   = Column(String(17), primary_key=True)
    release_id = Column(String(4),  primary_key=True)
    parents    = Column(Text)

    def __repr__(self):
        return "<NR_parents('%s','%s','%s')>" % (self.motif_id, self.release_id, self.parents)


class NR_release_diff(Base):
    """
    """
    __tablename__ = 'nr_release_diff'

    nr_release_id1 = Column(String(4), primary_key=True)
    nr_release_id2 = Column(String(4), primary_key=True)
    resolution     = Column(String(4), primary_key=True)
    direct_parent  = Column(Boolean)
    added_groups   = Column(Text)
    removed_groups = Column(Text)
    updated_groups = Column(Text)
    same_groups    = Column(Text)
    added_pdbs     = Column(Text)
    removed_pdbs   = Column(Text)
    num_added_groups   = Column(Integer)
    num_removed_groups = Column(Integer)
    num_updated_groups = Column(Integer)
    num_same_groups    = Column(Integer)
    num_added_pdbs     = Column(Integer)
    num_removed_pdbs   = Column(Integer)

    def __repr__(self):
        return "<NRReleaseDiff('%s','%s')>" % (self.nr_release_id1, self.nr_release_id2)


################################################################################
# Collections comparison classes
################################################################################
class NR_eqclass_collection:
    """
    A collection that is created either from a database release, or from a file.
    release = latest/previous/<empty>
    """

    def __init__(self, release='', file='', resolution=''):
        self.loops   = []
        self.groups  = []
        self.sl      = set()
        self.sg      = set()
        self.d       = collections.defaultdict(list) # [group: {loopX..loopY} ]
        self.sd      = collections.defaultdict(set)  # [group: set({loopX..loopY}]
        self.res     = resolution
        self.reps    = dict()
        self.release = release

        if release == '':
            pass
        elif release == 'latest':
            self.get_latest_release()
        elif release == 'previous':
            self.get_previous_release()
        else:
            self.release = release
            self.get_release()

        if file is not '':
            self.file = file
            self.read_motifs_csv_file()

        self.make_dictionary()
        self.make_sets()

    def read_motifs_csv_file(self):
        r = csv.reader(open(self.file), delimiter=',', quotechar='"')
        current_group = ''
        for row in r:
            self.loops.append(row[0])
            self.groups.append(row[1])
            if row[1] != current_group:
                current_group = row[1]
                self.reps[row[0]] = row[1]

    def make_sets(self):
        self.sl = set(self.loops)
        self.sg = set(self.groups)
        for k,v in self.d.iteritems():
            self.sd[k] = set(v)

    def make_dictionary(self):
        for i,loop in enumerate(self.loops):
            self.d[self.groups[i]].append(loop)

    def get_release(self):
        for loop in session.query(NR_pdb).filter(NR_pdb.release_id==self.release) \
                                         .filter(NR_pdb.class_id.like('NR_'+self.res+'_%')) \
                                         .all():
            self.loops.append(loop.id)
            self.groups.append(loop.class_id)

    def get_previous_release(self):
        if session.query(NR_release).first() is None:
            print 'No previous releases found'
            return
        release = session.query(NR_release).order_by(desc(NR_release.date))[0:2]
        if len(release) == 2:
            self.release = release[1].id
            self.get_release()
#             loops = session.query(NR_pdb).filter(NR_pdb.release_id==release[1].id) \
#                                          .filter(NR_pdb.class_id.like('%'+self.res+'%'))
#             for loop in loops:
#                 self.loops.append(loop.id)
#                 self.groups.append(loop.class_id)
        else:
            self.get_latest_release()


    def get_latest_release(self):
        if session.query(NR_release).first() is None:
            print 'No previous releases found'
            return
        release = session.query(NR_release).order_by(desc(NR_release.date)).first()
        self.release = release.id
        self.get_release()
#         loops = session.query(NR_pdb).filter(NR_pdb.release_id==release.id) \
#                                      .filter(NR_pdb.class_id.like('%'+self.res+'%'))
#         for loop in loops:
#             self.loops.append(loop.id)
#             self.groups.append(loop.class_id)


class NRCollectionMerger:
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
    When in release_diff mode, the class will not store the class ids, just the
    difference between releases
    release_mode: major/minor/reuse
    """
    def __init__(self, ensembles=None, release_mode='', release_description='', upload_mode=''):
        self.c            = ensembles # collections, NRCollectionMerger
        self.motifs       = []
        self.loops        = []
        self.history      = []
        self.final_ids    = dict()
        self.intersection = []
        self.release_diff = []
        self.added_groups       = []
        self.removed_groups     = []
        self.updated_groups     = []
        self.old_updated_groups = []
        self.same_groups        = []
        self.added_pdbs         = []
        self.removed_pdbs       = []

        self.upload_mode = upload_mode
        if upload_mode != 'release_diff':
            self.release = NR_release(mode=release_mode, description=release_description)
            self.__finalize()
            self.__process_set_diff()
            self.direct_parent = 1
            self.__process_release_diff()
            self.__commit()
        else:
            self.release = NR_release()
            self.release.id = self.c.c1.release
            self.direct_parent = 0
            self.__finalize()
            self.__process_release_diff()
            self.__release_diff_commit()


    def __finalize(self):
        """
        """
        self.added_pdbs   = self.c.c1.sl - self.c.c2.sl
        self.removed_pdbs = self.c.c2.sl - self.c.c1.sl

        for group_id in self.c.c1.sg:
            if group_id in self.c.new_ids:
                motif = NR_class(release_id=self.release.id,
                                 resolution=self.c.c1.res,
                                 comment=self.c.explanation[group_id])
                if self.c.parents.has_key(group_id):
                    parents = ','.join(set(self.c.parents[group_id]))
                else:
                    parents = ''
                self.added_groups.append(motif.id)
                print 'Group %s assigned new id %s' % (group_id, motif.id)

            elif group_id in self.c.correspond:
                old_id  = self.c.correspond[group_id]
                motif   = NR_class(id=old_id,
                                   release_id=self.release.id,
                                   increment=True,
                                   resolution=self.c.c1.res,
                                   comment=self.c.explanation[group_id])
                parents = ','.join(set([old_id] + self.c.parents[group_id]))
                self.updated_groups.append(motif.id)
                self.old_updated_groups.append(old_id)
                print 'Group %s corresponds to motif %s and is assigned new id %s' % (group_id, old_id, motif.id)

            elif group_id in self.c.exact_match:
                id = self.c.exact_match[group_id]
                motif = NR_class(id=id,
                                 release_id=self.release.id,
                                 resolution=self.c.c1.res,
                                 comment=self.c.explanation[group_id])
                parents = ''
                self.same_groups.append(motif.id)
                print 'Group %s matches exactly motif %s' % (group_id, motif.id)

            else:
                print 'Major problem'

            self.motifs.append(motif)
            self.final_ids[group_id] = motif.id

            if parents != '':
                self.history.append(NR_parents(class_id=motif.id, release_id=self.release.id, parents=parents))

            for loop_id in self.c.c1.d[group_id]:
                if loop_id in self.c.c1.reps:
                    self.loops.append(NR_pdb(id=loop_id, class_id=motif.id, release_id=self.release.id, rep = True))
                else:
                    self.loops.append(NR_pdb(id=loop_id, class_id=motif.id, release_id=self.release.id, rep = False))

            self.removed_groups = set(self.c.c2.groups) - set(self.old_updated_groups) - set(self.same_groups)


    def __process_release_diff(self):
        """
        """
        # do nothing during the very first upload
        if self.release.id == '0.1':
            return

        self.release_diff.append( NR_release_diff(
            nr_release_id1     = self.release.id,
            nr_release_id2     = self.c.c2.release,
            resolution         = self.c.c1.res,
            direct_parent      = self.direct_parent,
            added_groups       = ', '.join(self.added_groups),
            removed_groups     = ', '.join(self.removed_groups),
            updated_groups     = ', '.join(self.updated_groups),
            same_groups        = ', '.join(self.same_groups),
            added_pdbs         = ', '.join(self.added_pdbs),
            removed_pdbs       = ', '.join(self.removed_pdbs),
            num_added_groups   = len(self.added_groups),
            num_removed_groups = len(self.removed_groups),
            num_updated_groups = len(self.updated_groups),
            num_same_groups    = len(self.same_groups),
            num_added_pdbs     = len(self.added_pdbs),
            num_removed_pdbs   = len(self.removed_pdbs)
        ))


    def __process_set_diff(self):
        """
        """
        for loop_id in self.c.c1.sg:
            for motif_id in self.c.c2.sg:
                if motif_id != self.final_ids[loop_id] and \
                   self.c.intersection.has_key(loop_id) and \
                   self.c.intersection[loop_id].has_key(motif_id):

                    self.intersection.append( NR_setdiff(
                                nr_class1 = self.final_ids[loop_id],
                                nr_class2 = motif_id,
                                release_id = self.release.id,
                                intersection = ','.join(self.c.intersection[loop_id][motif_id]),
                                overlap = self.c.overlap[loop_id][motif_id],
                                one_minus_two = ','.join(self.c.setdiff[loop_id][motif_id]),
                                two_minus_one = ','.join(self.c.setdiff[motif_id][loop_id])
                    ))
                    self.intersection.append( NR_setdiff(
                                nr_class1 = motif_id,
                                nr_class2 = self.final_ids[loop_id],
                                release_id = self.release.id,
                                intersection = ','.join(self.c.intersection[loop_id][motif_id]),
                                overlap = self.c.overlap[motif_id][loop_id],
                                one_minus_two = ','.join(self.c.setdiff[motif_id][loop_id]),
                                two_minus_one = ','.join(self.c.setdiff[loop_id][motif_id])
                    ))


    def remove_release(self, release):
        session.query(NR_release).filter(NR_release.id==release).delete()
        session.query(NR_class).filter(NR_class.release_id==release).delete()
        session.query(NR_pdb).filter(NR_pdb.release_id==release).delete()
        session.query(NR_setdiff).filter(NR_setdiff.release_id==release).delete()
        session.query(NR_parents).filter(NR_parents.release_id==release).delete()
        session.query(NR_release_diff).filter(NR_release_diff.nr_release_id1==release).delete()
        session.query(NR_release_diff).filter(NR_release_diff.nr_release_id2==release).delete()

    def __release_diff_commit(self):
        """
        """
        try:
            session.add_all(self.release_diff)
            session.query(NR_handle).delete()
            session.commit()
            print 'Successful update'
        except sqlalchemy.exc.SQLAlchemyError, e:
            print 'Update failed. SQLAlchemy error. Rolling back.'
            print str(e)
            session.rollback()
            sys.exit()


    def __commit(self):
        """
        """
        try:
            handles = [i.handle for i in self.motifs]
            if len(handles) != len(set(handles)):
                pdb.set_trace()

            r = session.query(NR_release).filter(NR_release.id == self.release.id).first()
            if not r:
                session.add(self.release)

            session.add_all(self.motifs)
            session.add_all(self.loops)
            session.add_all(self.history)
            session.add_all(self.intersection)
            session.add_all(self.release_diff)

            session.commit()
            print 'Successful update'
        except sqlalchemy.exc.SQLAlchemyError, e:
            print 'Update failed. SQLAlchemy error. Rolling back.'
            print str(e)
            session.rollback()
            self.remove_release(self.release.id)
            sys.exit()
        except sqlalchemy.exc.DBAPIError, e:
            print 'Update failed. DBAPI error. Rolling back.'
            print str(e)
            session.rollback()
            self.remove_release(self.release.id)
            sys.exit()
        except sys.exc_info()[0]:
            print 'Update failed. Rolling back.'
            print sys.exc_info()[0]
            session.rollback()
            self.remove_release(self.release.id)
            sys.exit()

################################################################################
# End of database operations
################################################################################



Base.metadata.create_all(engine)
