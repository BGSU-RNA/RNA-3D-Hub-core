import logging
import collections
import csv

from sqlalchemy import desc

from pymotifs import core
from pymotifs.models import MlLoops
from pymotifs.models import MlReleases

logger = logging.getLogger(__name__)


class MotifCollection:
    def __init__(self, session_maker, file='', release='', type=''):
        self.loops = []
        self.groups = []
        self.type = type
        self.sl = set()
        self.sg = set()
        self.d = collections.defaultdict(list)  # [group: {loopX..loopY} ]
        self.sd = collections.defaultdict(set)  # [group: set({loopX..loopY}]
        self.release = release
        self.session = core.Session(session_maker)

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
        with open(self.file, 'rb') as raw:
            reader = csv.reader(raw, delimiter=',', quotechar='"')
            for row in reader:
                self.loops.append(row[0])
                self.groups.append(row[1])

    def make_sets(self):
        self.sl = set(self.loops)
        self.sg = set(self.groups)
        for k, v in self.d.iteritems():
            self.sd[k] = set(v)

    def make_dictionary(self):
        for i, loop in enumerate(self.loops):
            self.d[self.groups[i]].append(loop)

    def get_latest_release(self):
        with self.session() as session:
            current = session.query(MlReleases).\
                filter(MlReleases.type == self.type).\
                order_by(desc(MlReleases.date)).\
                first()

            if current is None:
                logger.info('No previous releases found')
                return

            self.release = current.id
        self.get_release()

    def get_release(self):
        with self.session() as session:
            query = session.query(MlLoops).\
                filter(MlLoops.release_id == self.release).\
                filter(MlLoops.type == self.type.upper())

            for loop in query:
                self.loops.append(loop.id)
                self.groups.append(loop.motif_id)
