"""

"""

import random, datetime, math, sys, pdb

from sqlalchemy import *
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker


engine  = create_engine('mysql://root:bioinfo@localhost/MotifVersions')
Session = sessionmaker(bind=engine)
session = Session()

Base = declarative_base()


def drop_all():
    Base.metadata.drop_all(engine) 


class Release(Base):
    """
    """
    __tablename__ = 'releases'
    
    id = Column(String(4), primary_key=True)
    date = Column(DateTime)
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
    __tablename__ = 'history'
    
    motif_id   = Column(String(11), primary_key=True)
    release_id = Column(String(4), 
                        ForeignKey('releases.id'),
                        primary_key=True )
    parents = Column(Text)  
    
    def __repr__(self):
        return "<Parents('%s','%s','%s')>" % (self.motif_id, self.release_id, self.parents)
    
#   release = relationship("Release", backref=backref('addresses'))
    
    
class Loop(Base):
    __tablename__ = 'loops'
    
    id         = Column(String(11), primary_key=True)
    motif_id   = Column(String(11), primary_key=True)
    release_id = Column(String(4), 
                        ForeignKey('releases.id'),
                        primary_key=True)

    def __repr__(self):
        return "<Loop('%s','%s','%s')>" % (self.id, self.motif_id, self.release_id)


class Motif(Base):
    """
    """
    __tablename__ = 'motifs'
    
    id         = Column(String(11), primary_key=True)
    release_id = Column(String(4), 
                        ForeignKey('releases.id'),
                        primary_key=True)
    type    = Column(String(2)) # IL, HL, JL
    handle  = Column(String(5)) # XXXXX
    version = Column(Integer)

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
    __tablename__ = 'loop_annotations'
    
    id          = Column(String(11), 
                         ForeignKey('loops.id'), 
                         primary_key=True)
    common_name = Column(Text)
    reference   = Column(Text)


class SetDiff(Base):
    __tablename__ = 'set_diff'
    
    motif_id1  = Column(String(11), primary_key=True)
    motif_id2  = Column(String(11), primary_key=True)
    release_id = Column(String(4),  primary_key=True)
    intersection  = Column(Text)
    overlap       = Column(String(4))
    one_minus_two = Column(Text)
    two_minus_one = Column(Text)

    def __repr__(self):
        return "<SetDiff('%s','%s','%s')>" % (self.motif_id1, self.motif_id2, self.release_id)        


class LoopCollection:
    """
    """
    def __init__(self, loops=[], groups=[], release=''):
        self.loops = loops
        self.groups = groups
        if release == '':
            pass
        elif release == 'latest':
            self.get_latest_release()
        else:
            self.release_id = release
            self.get_release()
        
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


class Uploader:
    """
    """
    def __init__(self, collections=None, mode='', description=''):
        self.c         = collections
        self.release   = Release(mode=mode, description=description)
        self.motifs    = []
        self.loops     = []
        self.history   = []
        self.intersection = []
        self.final_ids = dict()
    
        self.__finalize()
        self.__process_set_diff()
        self.show_report()
        self.__commit()

            
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
            elif group_id in self.c.correspond:
                old_id  = self.c.correspond[group_id]
                motif   = Motif(id=old_id, release_id=self.release.id, increment=True)
                parents = ','.join([old_id] + self.c.parents[group_id])
            elif group_id in self.c.exact_match:
                id = self.c.exact_match[group_id]
                motif = Motif(id=id, release_id=self.release.id)
                parents = ''
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
    
    def __commit(self):    
        try:
            session.add(self.release)
            session.add_all(self.motifs)
            session.add_all(self.loops)
            session.add_all(self.history)
            session.add_all(self.intersection)
            session.commit()
        except sys.exc_info()[0]:
            print 'Update failed. Rolling back.'
            print sys.exc_info()[0]
            session.rollback()
            self.remove_release(self.release.id)            
            
Base.metadata.create_all(engine)
            