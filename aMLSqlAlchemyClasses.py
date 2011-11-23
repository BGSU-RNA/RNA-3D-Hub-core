"""

"""

import random, datetime, math, sys, pdb, csv, os, shutil

from sqlalchemy import *
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

import sqlalchemy.exc

engine  = create_engine('mysql://root:bioinfo@localhost/MotifVersions')
Session = sessionmaker(bind=engine)
session = Session()

Base = declarative_base()


def drop_all():
    Base.metadata.drop_all(engine) 


class AllLoops(Base):
    """
    """
    __tablename__ = 'all_loops'
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
    loopName      = Column(Text)

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
        self.loopName  = parts[11]
                
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
    release_id     = Column(String(4))
        
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
    overlap       = Column(Float)
    one_minus_two = Column(Text)
    two_minus_one = Column(Text)

    def __repr__(self):
        return "<SetDiff('%s','%s','%s')>" % (self.motif_id1, self.motif_id2, self.release_id)        


class LoopOrder(Base):
    __tablename__ = 'loop_order'
        
    motif_id = Column(String(11), primary_key=True)
    loop_id  = Column(String(11), primary_key=True)
    release_id = Column(String(4), primary_key=True)
    original_order = Column(Integer)
    similarity_order = Column(Integer)

    def __repr(self):
        return "<LoopOrder('%s','%s','%i')>" % (self.motif_id, self.loop_id, self.original_order)


class LoopPosition(Base):
    __tablename__ = 'loop_positions'
    
    motif_id   = Column(String(11), primary_key=True)
    loop_id    = Column(String(11), primary_key=True)
    release_id = Column(String(4),  primary_key=True)
    nt_id      = Column(String(30), primary_key=True)
    position   = Column(Integer)
    
    def __repr__(self):
        return "<LoopPosition('%s','%s','%i')>" % (self.loop_id, self.nt_id, self.position)


class LoopDiscrepancy(Base):
    __tablename__ = 'mutual_discrepancy'
    
    loop_id1 = Column(String(11), primary_key=True)
    loop_id2 = Column(String(11), primary_key=True)
    release_id = Column(String(4), primary_key=True)
    discrepancy = Column(Float)
    
    def __repr__(self):
        return "<LoopDiscrepancy('%s','%s','%f')>" % (self.loop_id1, self.loop_id2, self.discrepancy)


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
    def __init__(self, collections=None, mode='', description='',files={}):
        self.c         = collections
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

        
        
Base.metadata.create_all(engine)
            