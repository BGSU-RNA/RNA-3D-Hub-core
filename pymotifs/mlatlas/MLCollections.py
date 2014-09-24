"""




"""

import logging
import collections
import csv

from sqlalchemy import desc


from models import Release, Loop, session

logger = logging.getLogger(__name__)


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
            logger.info('No previous releases found')
            return
        release = session.query(Release).\
                          filter(Release.type==self.type).\
                          order_by(desc(Release.date)).\
                          first()
        self.release = release.id
        self.get_release()

    def get_release(self):
        for loop in session.query(Loop).\
                            filter(Loop.release_id==self.release).\
                            filter(Loop.id.like(self.type+'%')).\
                            all():
            self.loops.append(loop.id)
            self.groups.append(loop.motif_id)
