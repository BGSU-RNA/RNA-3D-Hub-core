import sys
import csv
import random
import pdb
import MySQLdb as mdb

DATABASE = 'MotifAtlas'
HOST	 = 'localhost'
USER	 = 'root'
PASSWORD = 'rumata7'
FILE	 = 'input/motifs.csv'
MOTIFS_TABLE  = 'Motifs'
HISTORY_TABLE = 'History'


def connect_to_mysql():
	""" Establish database connection.
	"""
	try:
		conn = mdb.connect(HOST, USER, PASSWORD);
		cursor = conn.cursor()
		return cursor, conn
	except mdb.Error, e:  
		print "Unable to connect.\nError %d: %s" % (e.args[0],e.args[1])
		sys.exit(1)



def select_database(cursor):
	""" Use database. Create it if necessary.
	"""
	try:	
		cursor.execute("USE %s;" % DATABASE)
	except mdb.Error, e:
		print "No database %s. Trying to create it" % DATABASE
		cursor.execute("CREATE DATABASE %s;" % DATABASE)
		cursor.execute("USE %s;" % DATABASE)



def create_motifs_table(cursor):
	""" Make sure that the motifs table exists
	"""
	try:	
		cursor.execute("SELECT * FROM %s" % MOTIFS_TABLE)
	except mdb.Error,e:
		print "Table doesn't exist.\nError %d: %s" % (e.args[0],e.args[1])
		cursor.execute("""CREATE TABLE %s (loop_id VARCHAR(11),\
                          motif_id VARCHAR(11),version INT);""" % MOTIFS_TABLE)



def create_history_table(cursor):
	""" Make sure that the history table exists
	"""
	try:	
		cursor.execute("SELECT * FROM %s" % HISTORY_TABLE)
	except mdb.Error,e:
		print "Table doesn't exist.\nError %d: %s" % (e.args[0],e.args[1])
		cursor.execute("""CREATE TABLE %s (motif_id VARCHAR(11),\
                          parents TEXT,version INT);""" % HISTORY_TABLE)



def close_mysql_connection(cursor, conn):
	cursor.close()
	conn.close()



def read_motifs_csv_file():
	"""
	"""
	motifReader = csv.reader(open(FILE), delimiter=',', quotechar='"')
	loops = []
	groups = []
	for row in motifReader:
		loops.append(row[0])
		groups.append(row[1])
	return loops, groups



def get_version(cursor):
	""" Increment maximum existing version by 1. Set to 1 on the first run.
	"""	
	cursor.execute("SELECT MAX(version) FROM %s;" % MOTIFS_TABLE);
	oldVersion = cursor.fetchone()
	if oldVersion[0] is None:
		newVersion = 1
	else:
		newVersion = oldVersion[0] + 1
	print 'New version is %s' % newVersion
	return newVersion, oldVersion[0]



def get_previous_release(cursor,oldVersion):
	"""
	"""		
	cursor.execute("SELECT loop_id, motif_id FROM %s WHERE version = %i" 
	                % (MOTIFS_TABLE, oldVersion) )
	rows = cursor.fetchall()
	L = []; G = []
	for row in rows:
		L.append(row[0])
		G.append(row[1])
				
	return L, G				



def find_indices(tlist,m): 
	"""
	"""
	return [i for i,x in enumerate(tlist) if x == m]
				


def transform_data(loops, groups):
	"""
	"""	
	return dict([(i, [loops[x] for x in find_indices(groups,i)]) \
				for i in set(groups)])



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
	


def load_motifs_data(cursor, conn, M, version):
	""" 
	"""		
	for motif, loops in M.iteritems():
		for loop in loops:
			try:
				cursor.execute(("insert into %s" % MOTIFS_TABLE) + 
							   (" VALUES ('%s','%s',%s);" 
							   % (loop,motif,version)))
			except mdb.Error,e:
				print "Error %d: %s" % (e.args[0],e.args[1])
				sys.exit(1)				
	
	conn.commit()	
	print 'Motifs data loaded'



def load_history_data(cursor, conn, M, version):
	""" 
	"""		
	for child, parents in M.iteritems():
		concat = reduce(lambda x, y: x+y, parents)
		try:
			cursor.execute(("insert into %s" % HISTORY_TABLE) + 
						   (" VALUES ('%s','%s',%s);" 
						   % (child,concat,version)))
		except mdb.Error,e:
			print "Error %d: %s" % (e.args[0],e.args[1])
			sys.exit(1)
	
	conn.commit()	
	print 'History data loaded'



# def load_data_array(cursor, conn, loops, groups, version):
# 	""" 
# 	"""	
# 	M = dict([(x,generate_id(cursor,loops)) for x in set(groups)])
# 	
# 	for i,motif in enumerate(groups):
# 		try:
# 			cursor.execute(("insert into %s" % TABLE) + 
# 						   (" VALUES ('%s','%s',%s);" 
# 						   % (loops[i],M[groups[i]],version)))
# 		except mdb.Error,e:
# 			print "Error %d: %s" % (e.args[0],e.args[1])
# 	
# 	conn.commit()
# 	
# 	return M



def generate_id(cursor, loops):
	"""
	"""
	if loops[0][0:2] == 'IL':
		prefix = 'IL'
	elif loops[0][0:2] == 'HL':
		prefix = 'HL'
	else:
		print "%s Critical error: unable to identify id prefix." % loops[0][0:2]
		sys.exit(1)

	while True:
		id = '%s%05d' % (prefix, random.randint(1, 99999) )
		cursor.execute("SELECT motif_id FROM %s WHERE motif_id='%s';" % (MOTIFS_TABLE,id) )
		result = cursor.fetchone()
		if result is None:
			break
	
	return id


























