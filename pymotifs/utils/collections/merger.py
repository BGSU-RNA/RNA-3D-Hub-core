"""

The main logic of what groups are assigned new ids, are updated, or kept as is.

Analyzes two collections and establishes correspondences.
c1 = Collection1, new collection
c2 = Collection2, old collection

"""

import collections
import logging

logger = logging.getLogger(__name__)


class CollectionsMerger:

    def __init__(self, c1, c2, verbose=True):
        self.c1 = c1
        self.c2 = c2
        self.verbose = verbose

        self.intersection = collections.defaultdict(dict)
        self.noIntersection = dict()
        self.overlap = collections.defaultdict(dict)
        self.setdiff = collections.defaultdict(dict)

        self.correspond = dict()
        self.exact_match = dict()
        self.new_ids = []
        self.parents = collections.defaultdict(list)
        self.explanation = dict()

        self.minOverlap = float(2)/3

        self.compare_releases()
        self.show_report()
        self.establish_correspondences()

    def compare_releases(self):
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
                       if self.verbose:
                           logger.info('Groups %s and motif %s match up exactly' % (groupId, motifId))
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
                    for motifId, v in match.iteritems():
                        self.parents[groupId].append(motifId)
                        # if the overlap is big enough - assign the same id
                        if self.overlap[groupId][motifId] >= self.minOverlap and \
                           self.overlap[motifId][groupId] >= self.minOverlap:
                            success = True
                            self.correspond[groupId] = motifId
                            self.explanation[groupId] = 'Updated, 2 parents'
                            break
                    # if the overlap is not big enough - assign new id
                    if success is False:
                        self.new_ids.append(groupId)
                        self.explanation[groupId] = 'New id, 2 parents'
                    else:
                        # don't need to do anything because success==True
                        pass

                # with more than 2 motifs - assign new id
                elif len(match) > 2:
                    if self.verbose:
                        logger.info('Group %s has more than 2 parents',
                                    groupId)

                    self.explanation[groupId] = '> 2 parents'
                    self.new_ids.append(groupId)

                    for motifId, v in match.iteritems():
                        self.parents[groupId].append(motifId)

            # if didn't intersect with anything
            else:
                self.explanation[groupId] = 'New id, no parents'
                self.new_ids.append(groupId)

    def show_report(self):
        logger.info('%i new loops', len(self.c1.loops))
        logger.info('%i old loops', len(self.c2.loops))
        logger.info('%i new groups', len(self.c1.d.keys()))
        logger.info('%i old groups', len(self.c2.d.keys()))

        logger.info('%i new loops not in the database',
                    len(self.c1.sl - self.c2.sl))
        logger.info('%i old loops not in the new list',
                    len(self.c2.sl - self.c1.sl))

        all_new = set(self.c1.groups) - set(self.intersection.keys())
        if len(all_new) > 0:
            logger.info('%i entirely new groups', len(all_new))

        only_old = set(self.c2.groups) - set(self.intersection.keys())
        if len(only_old) > 0:
            logger.info('%i motifs only in the database', len(only_old))
