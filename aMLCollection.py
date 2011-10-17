"""

"""

import collections, csv, pdb


class MotifCollection:
    """
    
    """
    
    def __init__(self, file = '', collection = None):
        self.loops  = []
        self.groups = []
        self.sl     = set()
        self.sg     = set()        
        self.d      = collections.defaultdict(list) # [group: {loopX..loopY} ]
        self.sd     = collections.defaultdict(set)  # [group: set({loopX..loopY}]
        
        if collection is not None:
            self.loops  = collection.loops
            self.groups = collection.groups        
        elif file is not '':
            self.file = file
            self.read_motifs_csv_file()

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


class MotifCollectionMerger:
    """
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
                       print 'Groups %s and motif %s match up exactly' % (groupId, motifId)
                    # and the match is not perfect   
                    else:
                        self.parents[groupId].append(motifId)                    
                        if self.overlap[groupId][motifId] >= self.minOverlap and \
                           self.overlap[motifId][groupId] >= self.minOverlap:
                            self.correspond[groupId] = motifId
                            print 'Groups %s and motif %s match up with overlap' % (groupId, motifId)
                        else:
                            self.new_ids.append(groupId)
                            print 'Groups %s and motif %s have insufficient overlap' % (groupId, motifId)
                
                # with more than 2 motifs - assign new id
                elif len(match) > 2:
                    print 'Group %s has more than 2 parents' % groupId
                    self.new_ids.append(groupId)
                    for motifId,v in match.iteritems():
                        self.parents[groupId].append(motifId)
                
                # with one or two motifs
                else:
                    success = False
                    for motifId,v in match.iteritems():
                        self.parents[groupId].append(motifId)
                        # the overlap is big enough - assign the same id
                        if self.overlap[groupId][motifId] >= self.minOverlap and \
                           self.overlap[motifId][groupId] >= self.minOverlap:
                            success = True                              
                            self.correspond[groupId] = motifId
                            break
                    # if the overlap is not big enough - assign new id
                    if success == False:
                        self.new_ids.append(groupId)
                    else:
                        # don't need to do anything because success==True
                        pass
                        
            # if didn't intersect with anything
            else: 
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
        else: 
            print 'All groups matched some motifs in the database'      
              
        only_old = set(self.c2.groups) - set(self.intersection.keys())
        if len(only_old) > 0:
            print len(only_old), ' motifs only in the database'
            print [x for x in only_old]
        else: 
            print 'All motifs in the database matched some new groups'
            
