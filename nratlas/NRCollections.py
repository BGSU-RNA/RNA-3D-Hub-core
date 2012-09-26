"""

Classes for storing and comparing non-redundant lists.
Responsible for deciding what equivalence classes are new or should be updated.

"""

import logging
import collections
import csv

from sqlalchemy import desc


from NRSqlAlchemyClasses import NR_pdb, NR_release, session


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
                       logging.info('Groups %s and motif %s match up exactly' % (groupId, motifId))
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
                    logging.info('Group %s has more than 2 parents' % groupId)
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
        logging.info('%i new loops'  % len(self.c1.loops))
        logging.info('%i old loops'  % len(self.c2.loops))
        logging.info('%i new groups' % len(self.c1.d.keys()))
        logging.info('%i old groups' % len(self.c2.d.keys()))

        logging.info('%i new loops not in the database' % len(self.c1.sl - self.c2.sl))
        logging.info('%i old loops not in the new list' % len(self.c2.sl - self.c1.sl))

        all_new = set(self.c1.groups) - set(self.intersection.keys())
        if len(all_new) > 0:
            logging.info('%i entirely new groups' % len(all_new))
            print [x for x in all_new]

        only_old = set(self.c2.groups) - set(self.intersection.keys())
        if len(only_old) > 0:
            logging.info('%i motifs only in the database' % len(only_old))
            print [x for x in only_old]
