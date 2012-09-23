"""

"""

import random, datetime, math, sys, pdb, csv, os, shutil, collections, os
import logging

from sqlalchemy import *
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy.dialects.mysql import LONGTEXT, VARCHAR
from sqlalchemy.sql import or_

import sqlalchemy.exc
import ConfigParser

def get_engine():
    """
        looks up the connection parameters in a config file, which must be
        located in the same directory as the python scripts
    """
    script_path = os.path.abspath(os.path.dirname(sys.argv[0]))
    configfile = os.path.join(script_path, 'motifatlas.cfg')
    config = ConfigParser.RawConfigParser()
    config.read(configfile)

    return create_engine('mysql://' + config.get('database','user')     + ':' +
                                      config.get('database','password') + '@' +
                                      config.get('database','host')     + '/' +
                                      config.get('database','database'))
    return engine

engine = get_engine()

Session = sessionmaker(bind=engine)
session = Session()

Base = declarative_base()


################################################################################
# Motif tables declarations
################################################################################

class PdbBestChainsAndModels(Base):
    """
    """
    __tablename__  = 'pdb_best_chains_and_models'
    __table_args__ = ( UniqueConstraint('pdb_id'), )

    id          = Column(Integer, primary_key=True, autoincrement=True)
    pdb_id      = Column(String(4), index=True)
    best_chains = Column(VARCHAR(20, binary=True)) # case-sensitive
    best_models = Column(String(10))

class RedundantNucleotide(Base):
    """
    """
    __tablename__ = 'pdb_redundant_nucleotides'
    __table_args__ = ( UniqueConstraint('nt_id1', 'nt_id2'), )

    id     = Column(Integer, primary_key=True, autoincrement=True)
    pdb_id = Column(String(4), index=True)
    nt_id1 = Column(VARCHAR(30, binary=True))
    nt_id2 = Column(VARCHAR(30, binary=True))

class PdbObsolete(Base):
    """
    """
    __tablename__ = 'pdb_obsolete'
    __table_args__ = ( UniqueConstraint('obsolete_id'), )

    obsolete_id = Column(String(4), primary_key=True)
    replaced_by = Column(String(40), index=True)
    date = Column(DateTime)

class LoopPositions(Base):
    """
    """
    __tablename__ = 'loop_positions'
    __table_args__ = ( UniqueConstraint('loop_id', 'position'), )

    id       = Column(Integer, primary_key=True, autoincrement=True)
    loop_id  = Column(String(11))
    position = Column(Integer)
    nt_id    = Column(Text)
    bulge    = Column(Boolean)
    flanking = Column(Boolean)

class LoopSearchQA(Base):
    """
    """
    __tablename__ = 'loop_search_qa'
    __table_args__ = ( UniqueConstraint('loop_id1', 'loop_id2'), )

    id       = Column(Integer, primary_key=True, autoincrement=True)
    loop_id1 = Column(String(11))
    loop_id2 = Column(String(11))
    status   = Column(Integer)
    message  = Column(Text)


class LoopSearch(Base):
    """Stores information about pairwise all-against-all FR3D searches between
    all loops."""
    __tablename__ = 'loop_searches'
    __table_args__ = ( UniqueConstraint('loop_id1', 'loop_id2'), )

    id       = Column(Integer, primary_key=True, autoincrement=True)
    loop_id1 = Column(String(11))
    loop_id2 = Column(String(11))
    disc     = Column(Float)
    nt_list1 = Column(Text)
    nt_list2 = Column(Text)


class PairwiseInteractions(Base):
    """
    """
    __tablename__ = 'pdb_pairwise_interactions'

    iPdbSig    = Column(String(30), primary_key=True)
    jPdbSig    = Column(String(30), primary_key=True)
    pdb_id     = Column(String(4))
    f_lwbp     = Column(String(4))
    f_stacks   = Column(String(4))
    f_bphs     = Column(String(4))
    f_brbs     = Column(String(4))
    m_lwbp     = Column(String(3))
    m_mclw     = Column(String(10))
    m_mcOther  = Column(String(10))
    m_stacks   = Column(String(3))
    m_mcstacks = Column(String(8))
    r_sanger   = Column(String(6))
    r_lwbp     = Column(String(3))
    r_tertiary = Column(String(10))
    r_stacks   = Column(String(10))

    def __init__(self, iPdbSig='', jPdbSig=''):
        self.iPdbSig = iPdbSig
        self.jPdbSig = jPdbSig


class PdbInfo(Base):
    """
    """
    __tablename__ = 'pdb_info'
    __table_args__ = ( UniqueConstraint('structureId', 'chainId'), )

    id = Column(Integer, primary_key=True, autoincrement=True)
    structureId = Column(String(4), index=True)
    chainId = Column(VARCHAR(1, binary=True), index=True) # NB! case sensitive
    structureTitle = Column(Text)
    experimentalTechnique = Column(Text)
    depositionDate = Column(Date)
    releaseDate = Column(Date)
    revisionDate = Column(Text)
    ndbId = Column(Text)
    resolution = Column(Float)
    classification = Column(Text)
    structureMolecularWeight = Column(Float)
    macromoleculeType = Column(Text)
    structureAuthor = Column(Text)
    entityId = Column(Integer)
    sequence = Column(LONGTEXT)
    chainLength = Column(Integer)
    db_id = Column(Text)
    db_name = Column(Text)
    molecularWeight = Column(Float)
    secondaryStructure = Column(Text)
    entityMacromoleculeType = Column(Text)
#     ligandId = Column(Text)
#     ligandIdImage = Column(Text)
#     ligandMolecularWeight = Column(Float)
#     ligandFormula = Column(Text)
#     ligandName = Column(Text)
#     ligandSmiles = Column(Text)
#     InChI = Column(Text)
#     InChIKey = Column(Text)
    hetId = Column(Text)
    Ki = Column(Text)
    Kd = Column(Text)
    EC50 = Column(Text)
    IC50 = Column(Text)
    deltaG = Column(Text)
    deltaH = Column(Text)
    deltaS = Column(Text)
    Ka = Column(Text)
    compound = Column(Text)
    plasmid = Column(Text)
    source = Column(Text)
    taxonomyId = Column(Text)
    biologicalProcess = Column(Text)
    cellularComponent = Column(Text)
    molecularFunction = Column(Text)
    ecNo = Column(Text)
    expressionHost = Column(Text)
    cathId = Column(Text)
    cathDescription = Column(Text)
    scopId = Column(Text)
    scopDomain = Column(Text)
    scopFold = Column(Text)
    pfamAccession = Column(Text)
    pfamId = Column(Text)
    pfamDescription = Column(Text)
    crystallizationMethod = Column(Text)
    crystallizationTempK = Column(Float)
    phValue = Column(Float)
    densityMatthews = Column(Float)
    densityPercentSol = Column(Float)
    pdbxDetails = Column(Text)
    unitCellAngleAlpha = Column(Float)
    unitCellAngleBeta = Column(Float)
    unitCellAngleGamma = Column(Float)
    spaceGroup = Column(Text)
    lengthOfUnitCellLatticeA = Column(Float)
    lengthOfUnitCellLatticeB = Column(Float)
    lengthOfUnitCellLatticeC = Column(Float)
    Z_PDB = Column(Integer)
    rObserved = Column(Float)
    rAll = Column(Float)
    rWork = Column(Float)
    rFree = Column(Float)
    refinementResolution = Column(Float)
    highResolutionLimit = Column(Float)
    reflectionsForRefinement = Column(Integer)
    structureDeterminationMethod = Column(Text)
    conformerId = Column(Text)
    selectionCriteria = Column(Text)
#     fieldStrength = Column(Text)
#     manufacturer = Column(Text)
#     model = Column(Text)
    contents = Column(Text)
    solventSystem = Column(Text)
    ionicStrength = Column(Text)
    ph = Column(Text)
    pressure = Column(Text)
    pressureUnits = Column(Text)
    temperature = Column(Text)
    softwareAuthor = Column(Text)
    softwareName = Column(Text)
    version = Column(Text)
    method = Column(Text)
    details = Column(Text)
    conformerSelectionCriteria = Column(Text)
    totalConformersCalculated = Column(Text)
    totalConformersSubmitted = Column(Text)
    emdbId = Column(Text)
    emResolution = Column(Text)
    aggregationState = Column(Text)
    symmetryType = Column(Text)
    reconstructionMethod = Column(Text)
    specimenType = Column(Text)


class Dcc(Base):
    """
    """
    __tablename__ = 'dcc_residues'
    id                               = Column(String(20), primary_key=True)
    sfcheck_correlation              = Column(Float)
    sfcheck_correlation_side_chain   = Column(Float)
    sfcheck_real_space_R             = Column(Float)
    sfcheck_real_space_R_side_chain  = Column(Float)
    sfcheck_connect                  = Column(Float)
    sfcheck_shift                    = Column(Float)
    sfcheck_shift_side_chain         = Column(Float)
    sfcheck_density_index_main_chain = Column(Float)
    sfcheck_density_index_side_chain = Column(Float)
    sfcheck_B_iso_main_chain         = Column(Float)
    sfcheck_B_iso_side_chain         = Column(Float)
    mapman_correlation               = Column(Float)
    mapman_real_space_R              = Column(Float)
    mapman_Biso_mean                 = Column(Float)
    mapman_occupancy_mean            = Column(Float)


class Coordinates(Base):
    """
    """
    __tablename__ = 'pdb_coordinates'
    id          = Column(VARCHAR(30, binary=True), primary_key=True)
    pdb         = Column(String(4))
    pdb_type    = Column(String(4))
    model       = Column(Integer)
    chain       = Column(String(1))
    number      = Column(Integer)
    unit        = Column(String(3))
    ins_code    = Column(String(1))
    index       = Column(Integer)
    coordinates = Column(Text)


class Distances(Base):
    """
    """
    __tablename__ = 'pdb_distances'
    id1      = Column(VARCHAR(30, binary=True), primary_key=True, index=True)
    id2      = Column(VARCHAR(30, binary=True), primary_key=True, index=True)
    distance = Column(Float)


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

    def __repr__(self):
        return "<Loop('%s','%s')>" % (self.id, self.seq)


class PdbAnalysisStatus(Base):
    """Keeps track of what pdb files were analyzed and when. Will help to
       determine what programs need to be rerun if the algorithms or the
       annotations change. `pdb_file` is added to the primary key in order to
       keep track of asymmetric units and biological assemblies."""
    __tablename__ = 'pdb_analysis_status'
    id            = Column(String(4),  primary_key=True)
#     pdb_file    = Column(String(10), primary_key=True)
    distances     = Column(DateTime)
    coordinates   = Column(DateTime)
    interactions  = Column(DateTime)
    il            = Column(DateTime)
    hl            = Column(DateTime)
    j3            = Column(DateTime)
    qa            = Column(DateTime)
    motifs        = Column(DateTime)
    redundant_nts = Column(DateTime)
    best_chains_and_models = Column(DateTime)


class LoopQA(Base):
    """
    1 - valid loop
    2 - missing nucleotides
    3 - modified nucleotides
    4 - abnormal chain number
    5 - incomplete nucleotides
    6 - self-complementary internal loop
    """
    __tablename__ = 'loop_qa'
    id     = Column(String(11), primary_key=True)
    status = Column(Integer)
    modifications = Column(Text)
    nt_signature  = Column(Text)
    complementary = Column(Text)
    release_id = Column(String(4), primary_key=True)

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
    type        = Column(String(2), primary_key=True)
    date        = Column(DateTime)
    description = Column(Text)
    graphml     = Column(LONGTEXT)

    def __init__(self, mode='', description='', type=''):
        self.description = description
        self.mode = mode
        self.type = type
        self.graphml = ''
        self.compute_new_release_id()
        self.get_date()

    def __repr__(self):
        return "<Release('%s','%s','%s')>" % (self.id, self.date, self.description)

    def compute_new_release_id(self):
        prev = session.query(Release).filter(Release.type==self.type).order_by(desc(Release.date)).first()
        if prev is None:
            self.id = '0.1'
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

class ML_handle(Base):
    """
    """
    __tablename__ = 'ml_handles'
    id = Column(String(5), primary_key=True)


class Motif(Base):
    """
    """
    __tablename__ = 'ml_motifs'

    id         = Column(String(11), primary_key=True)
    release_id = Column(String(4),  primary_key=True)
    type       = Column(String(2)) # IL, HL, JL
    handle     = Column(String(5)) # XXXXX
    version    = Column(Integer)
    comment    = Column(Text)

    def __init__(self, id='', release_id='', increment=False, comment='', type=''):
        self.release_id = release_id
        self.comment    = comment
        self.type       = type
        if id == '':
            self.get_new_motif_id()
        elif increment is True:
            self.increment_motif_id(id)
        else:
            self.populate_fields(id)

    def __repr__(self):
        return "<Motif('%s','%s','%s')>" % (self.id, self.type, self.release_id)

    def get_new_motif_id(self):
#         self.type = 'IL'
        while True:
            self.handle = '%05d' % random.randrange(99999)
            if session.query(Motif).filter(Motif.handle==self.handle).first() is None:
                if session.query(ML_handle).filter(ML_handle.id==self.handle).first() is None:
                    h = ML_handle(id=self.handle)
                    session.add(h)
                    break

        self.version = 1
        self.id = self.type + '_' + self.handle + '.' + str(self.version)

    def increment_motif_id(self, id):
#         self.type = 'IL'
        self.handle  = id[3:8]
        self.version = int(id[9:]) + 1
        self.id = self.type + '_' + self.handle + '.' + str(self.version)

    def populate_fields(self, id):
        self.id = id
#         self.type = 'IL'
        self.handle  = id[3:8]
        self.version = int(id[9:])



class MotifAnnotation(Base):
    __tablename__ = 'ml_motif_annotations'

    motif_id     = Column(String(11), primary_key=True)
    common_name  = Column(Text)
    annotation   = Column(Text)
    author       = Column(String(20))
    bp_signature = Column(String(100))
    date         = Column(DateTime)


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

class Release_diff(Base):
    """
    """
    __tablename__ = 'ml_release_diff'

    release_id1    = Column(String(4), primary_key=True)
    release_id2    = Column(String(4), primary_key=True)
    type           = Column(String(2), primary_key=True)
    direct_parent  = Column(Boolean)
    added_groups   = Column(Text)
    removed_groups = Column(Text)
    updated_groups = Column(Text)
    same_groups    = Column(Text)
    added_loops    = Column(Text)
    removed_loops  = Column(Text)
    num_added_groups   = Column(Integer)
    num_removed_groups = Column(Integer)
    num_updated_groups = Column(Integer)
    num_same_groups    = Column(Integer)
    num_added_loops    = Column(Integer)
    num_removed_loops  = Column(Integer)

    def __repr__(self):
        return "<NRReleaseDiff('%s','%s')>" % (self.release_id1, self.release_id2)



################################################################################
# End of motif tables declarations
################################################################################

################################################################################
# Collection comparison classes
################################################################################

class MotifCollection:
    """

    """
    def __init__(self, file='', release='',type=''):
        self.loops   = []
        self.groups  = []
        self.type    = type
        self.sl      = set()
        self.sg      = set()
        self.d       = collections.defaultdict(list) # [group: {loopX..loopY} ]
        self.sd      = collections.defaultdict(set)  # [group: set({loopX..loopY}]
        self.release = release

        if file is not '':
            self.file = file
            self.read_motifs_csv_file()
        elif release == 'latest':
            self.get_latest_release()
        elif release != '':
            self.release = release
            self.get_release()
        else:
            pass

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
        if session.query(Release).filter(Release.type==self.type).first() is None:
            print 'No previous releases found'
            return
        release = session.query(Release).filter(Release.type==self.type).order_by(desc(Release.date)).first()
        self.release = release.id
        self.get_release()

    def get_release(self):
        for loop in session.query(Loop).filter(Loop.release_id==self.release).filter(Loop.id.like(self.type+'%')).all():
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
    def __init__(self, ensembles=None, mode='', description='', files={}, upload_mode='', motif_type=''):
        self.c           = ensembles
        self.motif_type  = motif_type.upper()
        self.release     = Release(mode=mode, description=description, type=self.motif_type)
        self.upload_mode = upload_mode
        self.motifs    = []
        self.loops     = []
        self.history   = []
        self.intersection = []
        self.final_ids    = dict()
        self.files        = files

        self.release_diff       = []
        self.added_groups       = []
        self.removed_groups     = []
        self.updated_groups     = []
        self.old_updated_groups = []
        self.same_groups        = []
        self.added_loops         = []
        self.removed_loops       = []

        self.loop_order       = []
        self.loop_positions   = []
        self.loop_discrepancy = []
        self.bp_signatures = []

    def import_release(self):
        """
        """
        if self.upload_mode != 'release_diff':
            self.direct_parent = 1
            self.__finalize()
            self.__process_set_diff()
            self.__process_release_diff()
            self.__process_motif_loop_order()
            self.__process_motif_loop_positions()
            self.__process_mutual_discrepancy()
            self.__process_bp_signatures()
            self.__store_correspondences()
#             self.show_report()
            self.__commit()
            self.__rename_mat_files()
            self.__move_2d_files()
            self.__process_graphml_file()
        else:
            self.direct_parent = 0
            self.release = Release()
            self.release.id = self.c.c1.release
            self.__finalize()
            self.__process_release_diff()
            self.__release_diff_commit()


    def __finalize(self):
        """
        """
        self.added_loops   = self.c.c1.sl - self.c.c2.sl
        self.removed_loops = self.c.c2.sl - self.c.c1.sl

        for group_id in self.c.c1.sg:
            if group_id in self.c.new_ids:
                motif = Motif(release_id=self.release.id,
                              comment=self.c.explanation[group_id],
                              type=self.motif_type)
                self.added_groups.append(motif.id)
                if self.c.parents.has_key(group_id):
                    parents = ','.join(self.c.parents[group_id])
                else:
                    parents = ''
                print 'Group %s assigned new id %s' % (group_id, motif.id)
            elif group_id in self.c.correspond:
                old_id  = self.c.correspond[group_id]
                motif   = Motif(id=old_id,
                                release_id=self.release.id,
                                increment=True,
                                comment=self.c.explanation[group_id],
                                type=self.motif_type)
                parents = ','.join(set([old_id] + self.c.parents[group_id]))
                self.updated_groups.append(motif.id)
                self.old_updated_groups.append(old_id)
                print 'Group %s corresponds to motif %s and is assigned new id %s' % (group_id, old_id, motif.id)
            elif group_id in self.c.exact_match:
                id = self.c.exact_match[group_id]
                motif = Motif(id=id,
                              release_id=self.release.id,
                              comment=self.c.explanation[group_id],
                              type=self.motif_type)
                self.same_groups.append(motif.id)
                parents = ''
                print 'Group %s matches exactly motif %s' % (group_id, motif.id)
            else:
                print 'Major problem'

            self.motifs.append(motif)
            if parents != '':
                if self.upload_mode != 'release_diff':
                    self.__inherit_history(motif, parents)
                self.history.append(Parents(motif_id=motif.id, release_id=self.release.id, parents=parents))
            self.final_ids[group_id] = motif.id
            for loop_id in self.c.c1.d[group_id]:
                self.loops.append(Loop(id=loop_id, motif_id=motif.id, release_id=self.release.id))

            self.removed_groups = set(self.c.c2.groups) - set(self.old_updated_groups) - set(self.same_groups)

    def __store_correspondences(self):
        """
        """
        f = open(self.files['correspondences'], 'w')
        for group_id, motif_id in self.final_ids.iteritems():
            f.write('%s,%s\n' % (group_id, motif_id) )
        f.close()

    def __inherit_history(self, motif, parents):
        """
        """
        parent_motifs = parents.split(',')
        common_name = []
        annotation = []
        author = []
        for parent in parent_motifs:
            parent_motif = session.query(MotifAnnotation).filter(MotifAnnotation.motif_id==parent).first()
            if parent_motif:
                if parent_motif.common_name is not None and parent_motif.common_name != '':
                    common_name.append(parent_motif.common_name + ' [' + parent + ']')
                if parent_motif.annotation is not None and parent_motif.annotation != '':
                    annotation.append(parent_motif.annotation + ' [' + parent + ']')
                if parent_motif.author is not None and parent_motif.author != '':
                    author.append(parent_motif.author)

        session.merge(MotifAnnotation(motif_id=motif.id,
                                      common_name=', '.join(common_name),
                                      annotation=', '.join(annotation),
                                      author=', '.join(author)))

    def __process_release_diff(self):
        """
        """
        # do nothing during the very first upload
        if self.release.id == '0.1':
            return

        self.release_diff.append( Release_diff(
            release_id1        = self.release.id,
            release_id2        = self.c.c2.release,
            type               = self.motif_type,
            direct_parent      = self.direct_parent,
            added_groups       = ', '.join(self.added_groups),
            removed_groups     = ', '.join(self.removed_groups),
            updated_groups     = ', '.join(self.updated_groups),
            same_groups        = ', '.join(self.same_groups),
            added_loops         = ', '.join(self.added_loops),
            removed_loops       = ', '.join(self.removed_loops),
            num_added_groups   = len(self.added_groups),
            num_removed_groups = len(self.removed_groups),
            num_updated_groups = len(self.updated_groups),
            num_same_groups    = len(self.same_groups),
            num_added_loops     = len(self.added_loops),
            num_removed_loops   = len(self.removed_loops)
        ))

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
        """
           Must pay attention to release type so that only IL or HL release
           is deleted, not both
        """
# DELETE FROM ml_set_diff WHERE release_id='0.7';
# DELETE FROM ml_mutual_discrepancy WHERE release_id='0.7';
# DELETE FROM ml_motifs WHERE release_id='0.7';
# DELETE FROM ml_loops WHERE release_id='0.7';
# DELETE FROM ml_loop_positions WHERE release_id='0.7';
# DELETE FROM ml_loop_order WHERE release_id='0.7';
# DELETE FROM ml_history WHERE release_id='0.7';
# DELETE FROM ml_releases WHERE id='0.7';
# DELETE FROM ml_release_diff WHERE release_id1='0.7' OR `release_id2`='0.7';

        # ml_handles
        logging.info('Removing release %s from table %s' % (release, ML_handle.__tablename__) )
        session.query(ML_handle).delete(synchronize_session='fetch')
        # ml_loop_order
        logging.info('Removing release %s from table %s' % (release, LoopOrder.__tablename__) )
        session.query(LoopOrder).\
                filter(LoopOrder.release_id==release).\
                filter(LoopOrder.motif_id.like(self.motif_type+'%')).\
                delete(synchronize_session='fetch')
        # ml_loop_positions
        logging.info('Removing release %s from table %s' % (release, LoopPosition.__tablename__) )
        session.query(LoopPosition).\
                filter(LoopPosition.release_id==release).\
                filter(LoopPosition.motif_id.like(self.motif_type+'%')).\
                delete(synchronize_session='fetch')
        # ml_releases
        logging.info('Removing release %s from table %s' % (release, Release.__tablename__) )
        session.query(Release).\
                filter_by(id=release).\
                filter_by(type=self.motif_type).\
                delete(synchronize_session='fetch')
        # ml_motifs
        logging.info('Removing release %s from table %s' % (release, Motif.__tablename__) )
        session.query(Motif).\
                filter_by(release_id=release).\
                filter_by(type=self.motif_type).\
                delete(synchronize_session='fetch')
        # ml_loops
        logging.info('Removing release %s from table %s' % (release, Loop.__tablename__) )
        session.query(Loop).\
                filter(Loop.release_id==release).\
                filter(Loop.id.like(self.motif_type+'%')).\
                delete(synchronize_session='fetch')
        # ml_set_diff
        logging.info('Removing release %s from table %s' % (release, SetDiff.__tablename__) )
        session.query(SetDiff).\
                filter(SetDiff.release_id==release).\
                filter(SetDiff.motif_id1.like(self.motif_type+'%')).\
                delete(synchronize_session='fetch')
        # ml_history
        logging.info('Removing release %s from table %s' % (release, Parents.__tablename__) )
        session.query(Parents).\
                filter(Parents.release_id==release).\
                filter(Parents.motif_id.like(self.motif_type+'%')).\
                delete(synchronize_session='fetch')
        # ml_mutual_discrepancy
        logging.info('Removing release %s from table %s' % (release, LoopDiscrepancy.__tablename__) )
        session.query(LoopDiscrepancy).\
                filter(LoopDiscrepancy.release_id==release).\
                delete(synchronize_session='fetch')
        # ml_release_diff
        logging.info('Removing release %s from table %s' % (release, Release_diff.__tablename__) )
        session.query(Release_diff).\
                filter(or_(Release_diff.release_id1==release, Release_diff.release_id2==release)).\
                filter_by(type=self.motif_type).\
                delete(synchronize_session='fetch')
        session.commit()

    def __release_diff_commit(self):
        """
        """
        try:
            session.add_all(self.release_diff)
            session.query(ML_handle).delete()
            session.commit()
            print 'Successful update'
        except sqlalchemy.exc.SQLAlchemyError, e:
            print 'Update failed. SQLAlchemy error. Rolling back.'
            print str(e)
            session.rollback()
            sys.exit()

    def __commit(self):
        try:
            session.add(self.release)
            session.add_all(self.motifs)
            session.add_all(self.loops)
            session.add_all(self.history)
            session.add_all(self.intersection)
            session.add_all(self.release_diff)

            session.add_all(self.loop_order)
            session.add_all(self.loop_positions)
            session.add_all(self.loop_discrepancy)
            session.commit()
            print 'Successful update'
        except sqlalchemy.exc.SQLAlchemyError, e:
            print 'Update failed. Rolling back.'
            print str(e)
            session.rollback()
            self.remove_release(self.release.id)
            sys.exit()
        except sqlalchemy.exc.DBAPIError, e:
            print 'Update failed. Rolling back.'
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


    def __process_bp_signatures(self):
        """
        "Group_001","06_cWW-cWW-cWW"
        """
        r = csv.reader(open(self.files['BpSignatures']), delimiter=',', quotechar='"')
        for row in r:
            session.merge(MotifAnnotation(motif_id=self.final_ids[row[0]],
                                          bp_signature=row[1]))
        session.commit()


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

        if not os.path.exists(self.files['MatFiles_dest']):
            os.mkdir(self.files['MatFiles_dest'])

        for file in self.c.c1.sg:
            src = os.path.join(self.files['MatFiles_origin'], file+'.mat')
            dst = os.path.join(self.files['MatFiles_dest'], self.final_ids[file] + '.mat')
            if os.path.exists(src):
                shutil.copyfile(src, dst)
            else:
                print "File %s wasn't found" % src


    def __move_2d_files(self):
        """
        """
        img_release_folder = os.path.join(self.files['2ds_destination'],self.motif_type + self.release.id)
        if not os.path.exists(img_release_folder):
            os.mkdir(img_release_folder)

        for file in self.c.c1.sg:
            src = os.path.join(self.files['2ds_origin'], file + '.png')
            dst = os.path.join(img_release_folder, self.final_ids[file] + '.png')
            if os.path.exists(src):
                shutil.copyfile(src, dst)
            else:
                print "File %s wasn't found" % src


    def __process_graphml_file(self):
        """
        """
        graphml = os.path.join(self.files['folder'],'Supergroups.graphml')
        if os.path.exists(graphml):
            f = open(graphml, 'r')
            contents = f.read()
            for group in self.c.c1.sg:
#                 if group in contents:
#                     # first try replacing group names as they were given
#                     contents = contents.replace(group, self.final_ids[group])
#                 else:
#                     # try to replace ids like Group_001 in case the matfiles were renamed to IL_001 etc
                parts = group.split('_')
                contents = contents.replace('Group_'+parts[1], self.final_ids[group])
            contents = contents.replace('\n','')
            self.release.graphml = contents
            session.add(self.release)
            session.commit()


################################################################################
# End of database operations
################################################################################


Base.metadata.create_all(engine)
