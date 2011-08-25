def update_2d_dict(D,a,b,c):
	"""
	"""
	if a in D:
		D[a][b] = c
	else:
		D[a] = {b:c}
	return D
	
		
		
def check_for_duplicates(L):
	if len(L) != len(set(L)):
		print 'There are duplicates in the dataset'
		
		

def compare_releases(cursor, m, M, l, L):
	"""
	m - new motifs
	M - old motifs
	l - list of all new loops
	L - list of all loops already in the database
	"""
	Map = dict()      # final result
	Parents = dict()  # history	
	
	I = dict() # intersection with previous release, keys - group ids
	R = dict() # reverse intersection, keys - motif ids
	N = dict() # no intersection with previous release
	O = dict() # % overlap
	D = dict() # set difference between overlapping groups
	
	for groupId, group in m.iteritems():		
		a_set = set(group)
		for motifId, motif in M.iteritems():
			b_set = set(motif)
			inters = a_set.intersection(b_set)
			if len(inters) > 0:			
				I = update_2d_dict(I, groupId, motifId, inters)
				R = update_2d_dict(R, motifId, groupId, inters)
				O = update_2d_dict(O, groupId, motifId, 
								  [float(len(inters))/len(group), 
								  float(len(inters))/len(motif)])
				if O[groupId][motifId] != [100,100]:
					D = update_2d_dict(D, groupId, motifId, 
					                   [a_set - b_set, b_set - a_set])
			
	#pdb.set_trace()	

	print 'Total %i loops in the database' % len(L)
	print 'Total %i loops in the new groups' % len(l)
	print 'Total %i motifs in the database' % len(M)
	print 'Total %i groups in the new version' % len(m)
	print '++++++++++++++++++++++++++++++++++++++++++++'
	
	appeared    = set(l) - set(L)
	if len(appeared) > 0:
		print 'New loops that are not in the database:', appeared
	else:
		print 'All loops are already present in the database'
		
	disappeared = set(L) - set(l)
	if len(disappeared) > 0:
		print 'Obsoleted loops:', disappeared
	else:
		print 'No loops were obsoleted'
	print '++++++++++++++++++++++++++++++++++++++++++++'

	# groupIds that didn't intersect with any motifId
	entirelyNewGroups = set(m.keys()) - set(I.keys())
	if len(entirelyNewGroups) > 0:
		print 'Entirely new groups:'
		for i in entirelyNewGroups: print i
	else:
		print 'All groups matched some motifs in the database'

	# motifIds that didn't intersect with any groupId
	disappearedMotifs = set(M.keys()) - set(R.keys())
	if len(disappearedMotifs) > 0:
		print 'Motif groups that didn\'t match anything:'
		for i in disappearedMotifs: print i
	else:
		print 'All motifs in the database matched some new groups'
	print '++++++++++++++++++++++++++++++++++++++++++++'
	
	# assign new id if:
	# 	more than 2 parents
	#	intersection < 2/3
	#   the size of the final cluster is more than 30% larger
	
	for groupId, matches in I.iteritems():
	
		# matched only 1 motif group and that motif group didn't match anything else
 		if len(matches) == 1 and O[groupId].values() == [[1,1]]:
 			Map[O[groupId].keys()[0]] = m[groupId]
 			print 'Identical match: %s - %s' % (groupId, O[groupId].keys()[0])
 		
 		# more than 2 parents, create new id
		elif len(matches) > 2:
			print 'Group %s has more than 2 parents' % groupId
			newId = generate_id(cursor,m[groupId])
			Map[newId] = m[groupId]
			Parents[newId] = matches.keys()
		else:
			# loop over matched motif groups, select the one with the largest match
			MINOVERLAP = float(2)/3
			reuse = ''
			for motifId, overlap in O[groupId].iteritems():					
				if overlap[0] >= MINOVERLAP and overlap[1] >= MINOVERLAP:
					print """Intersection between group %s and motif %s is (%f,%f).\
					Reuse the motif id""" % (groupId, motifId, overlap[0], overlap[1])
					reuse = motifId
					break
					
			# if no motif id can be reused, create a new one
			if reuse != '':
				Map[reuse] = m[groupId]
				Parents[reuse] = matches.keys()
			else:
				newId = generate_id(cursor,m[groupId])
				print """Created a new id %s for group %s""" % (newId, groupId)
				print 'Intersection', I[groupId]
				print 'Overlap', O[groupId]
				print 'Set difference', D[groupId]
				Map[newId] = m[groupId]
				Parents[newId] = I[groupId].keys()
	 			
	pdb.set_trace()	
	
	return Map, Parents
