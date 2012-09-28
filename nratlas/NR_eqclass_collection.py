"""

Class for storing and comparing non-redundant lists.

"""

import logging
import collections
import csv

from sqlalchemy import desc


from models import NR_pdb, NR_release, session


class NR_eqclass_collection:
    """
        A collection that is created either from a database release or a file.
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
            logging.info('No previous releases found')
            return
        release = session.query(NR_release).order_by(desc(NR_release.date))[0:2]
        if len(release) == 2:
            self.release = release[1].id
            self.get_release()
        else:
            self.get_latest_release()

    def get_latest_release(self):
        if session.query(NR_release).first() is None:
            logging.info('No previous releases found')
            return
        release = session.query(NR_release).order_by(desc(NR_release.date)).first()
        self.release = release.id
        self.get_release()
