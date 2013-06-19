"""

A script responsible for connecting to the database, managing the session,
creating all tables. Also contains declarative sqlalchemy table definitions.

"""


import random
import datetime
import math
import sys
import pdb
import os
import ConfigParser

import collections

from sqlalchemy import *
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy.dialects.mysql import LONGTEXT, VARCHAR
from sqlalchemy.sql import or_

import sqlalchemy.exc


Base = declarative_base()

# Utility functions
def drop_all():
    Base.metadata.drop_all(engine)

def create_all():
    Base.metadata.create_all(engine)

def get_session():
    """
        looks up the connection parameters in a config file, which must be
        located in the same directory as the python scripts.
    """

    filename = 'motifatlas.cfg'
    script_path = os.path.dirname(os.path.abspath( __file__ ))
    configfile = os.path.join(script_path, filename)
    config = ConfigParser.RawConfigParser()
    config.read(configfile)

    # use the test database when run from unittests
    if 'unittest' in sys.modules:
        env = 'test'
    else:
        env = config.get('general', 'environment')

    if env == 'test':
        db = 'database_test'
    else:
        db = 'database'

    print '%s environment' % env
    print 'Connecting to the `%s` database' % config.get('database', db)

    engine = create_engine('mysql://' + config.get('database','user')     + ':' +
                                        config.get('database','password') + '@' +
                                        config.get('database','host')     + '/' +
                                        config.get('database', db))

    Session = sessionmaker(bind=engine)
    return Session(), engine


session, engine = get_session()


class PdbUnitIdCorrespondence(Base):
    """
    """
    __tablename__  = 'pdb_unit_id_correspondence'
    __table_args__ = ( UniqueConstraint('unit_id', 'old_id'), )

    id       = Column(Integer, primary_key=True, autoincrement=True)
    old_id   = Column(VARCHAR(30, binary=True), index=True)
    unit_id  = Column(VARCHAR(30, binary=True), index=True)
    pdb      = Column(String(4), index=True)
    model    = Column(Integer)    # defaults to 1
    chain    = Column(VARCHAR(1, binary=True)) # chain id, case sensitive
    seq_id   = Column(Integer)    # residue number
    comp_id  = Column(String(3))  # component id: ALA, A, HOH
    atom     = Column(String(3))  # atom name. If blank, all atoms in residue
    alt_id   = Column(String(1))  # occupancy: blank, A, B, or 0
    ins_code = Column(String(1))  # A for residues like 109A
    sym_op   = Column(String(20)) # 1_555, 6_555
    pdb_file = Column(String(6))  # .pdb, .pdb10


class PdbBestChainsAndModels(Base):
    """
    """
    __tablename__  = 'pdb_best_chains_and_models'
    __table_args__ = ( UniqueConstraint('pdb_id'), )

    id          = Column(Integer, primary_key=True, autoincrement=True)
    pdb_id      = Column(String(4), index=True)
    best_chains = Column(VARCHAR(50, binary=True)) # case-sensitive
    best_models = Column(String(50))


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
    border   = Column(Boolean)


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
    f_bphs     = Column(String(5)) # n0BPh - can be 5 characters
    f_brbs     = Column(String(4))
    f_crossing = Column(Integer)
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
    unit_ids      = Column(DateTime)
    unit_ordering = Column(DateTime)


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
        return "<LoopQA('%s','%s','%s')>" % (self.id, self.status, self.release_id)


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
    annotation  = Column(Text)
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
    """equivalent to nr_parents"""
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
        self.handle  = id[3:8]
        self.version = int(id[9:]) + 1
        self.id = self.type + '_' + self.handle + '.' + str(self.version)

    def populate_fields(self, id):
        self.id = id
        self.handle  = id[3:8]
        self.version = int(id[9:])


class MotifAnnotation(Base):
    __tablename__ = 'ml_motif_annotations'

    motif_id     = Column(String(11), primary_key=True)
    common_name  = Column(Text)
    annotation   = Column(Text)
    author       = Column(Text)
    bp_signature = Column(Text)
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


class NR_pdb(Base): # = loop
    """
    """
    __tablename__ = 'nr_pdbs'
    id         = Column(String(4),  primary_key=True)
    class_id   = Column(String(17), primary_key=True)
    release_id = Column(String(6),  primary_key=True)
    rep        = Column(Boolean)

    def __repr__(self):
        return "<PDB('%s','%s','%s')>" % (self.id, self.class_id, self.release_id)


class NR_class(Base): # = motif
    """
    """
    __tablename__ = 'nr_classes'

    id         = Column(String(17), primary_key=True)
    release_id = Column(String(6),  primary_key=True)
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

    id          = Column(String(6), primary_key=True)
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
    release_id    = Column(String(6),  primary_key=True)
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
        return "<NR_parents('%s','%s','%s')>" % (self.class_id, self.release_id, self.parents)


class NR_release_diff(Base):
    """
    """
    __tablename__ = 'nr_release_diff'

    nr_release_id1 = Column(String(6), primary_key=True)
    nr_release_id2 = Column(String(6), primary_key=True)
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


class PdbUnitOrdering(Base):
    """
    Stores the ordering of units in a PDB file. The id used is an old style ID.
    """
    __tablename__ = 'pdb_unit_ordering'
    nt_id = Column(VARCHAR(30, binary=True), primary_key=True)
    pdb = Column(String(4), index=True)
    index = Column(Integer)

    def __repr__(self):
        return "<PdbUnitOrdering('%s','%s')>" % (self.nt_id, self.index)


class PdbModifiedCorrespondecies(Base):
    """
    Stores what modified bases correspond to which standard RNA bases.
    """
    __tablename__ = 'pdb_modified_correspondecies'
    id = Column(Integer, primary_key=True)
    modified_unit = Column(String(10))
    standard_unit = Column(String(1))


class Feature(Base):
    __tablename__ = 'feature_info'
    id = Column(Integer, primary_key=True)
    name = Column(String(50))
    pdb = Column(String(4))
    type_id = Column(Integer, ForeignKey('feature_type.id'))


class FeatureType(Base):
    __tablename__ = 'feature_type'
    id = Column(Integer, primary_key=True)
    name = Column(String(50))


class FeatureNucleotides(Base):
    __tablename__ = 'feature_nucleotides'
    id = Column(Integer, primary_key=True)
    unit_id = Column(VARCHAR(30, binary=True))
    feature_id = Column(Integer, ForeignKey('feature_info.id'))


class LoopLocations(Base):
    __tablename__ = 'loop_locations'
    name = Column(String(50), unique=True)
    featue_type_id = Column(Integer, ForeignKey('feature_type.id'))


class LoopLocationAnnotation(Base):
    __tablename__ = 'loop_location_annotation'
    loop_id = Column(String(11))
    location_id = Column(Integer, ForeignKey('loop_locations.id'))


Base.metadata.create_all(engine)
