"""

8/24/2011 Anton Petrov 
Dump the entire release history into a graphml file for visualization with Cytoscape

"""

import pdb
import aMLDbConnector as adb


FILE = 'graphml.txt'


def aMotifsToGraphML():

	db = adb.aDbConnector()	
	
	# get node info from the database
	rows = db.query("SELECT * FROM %s" % db.motifsTable)
	
	M = dict()
	V = dict()
	for row in rows:
		id = "%s.%i" % (row[1],row[2])
		if id not in M: M[id] = []
		M[id].append(row[0])
		if row[1] not in V: V[row[1]] = []
		V[row[1]].append(row[2])
	
	# get history info from the database
	rows = db.query("SELECT * FROM %s" % db.historyTable)
	
	P = dict()
	for row in rows:
		id = "%s.%i" % (row[0],row[2])
		t = ['IL'+x for x in row[1].split('IL')]
		t.remove('IL')
		P[id] = t
	
	# output graphml file
	f = open(FILE, 'w')
	f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
	f.write('<graphml xmlns="http://graphml.graphdrawing.org/xmlns">\n')
	f.write('<graph id="G" edgedefault="undirected">\n')
	f.write('<key id="invisible" for="node" attr.name="invisible" attr.type="double"/>\n')
	f.write('<key id="invisible" for="edge" attr.name="invisible" attr.type="double"/>\n')
	f.write("<node id='root'/>\n")	
	
	for motif in M:
		f.write("<node id='%s'/>\n" % motif)
	
	for motif in V:
		minVersion = min(V[motif])
		maxVersion = max(V[motif])
	
		if maxVersion > 1:
			if minVersion > 1:
				f.write("<edge source='root' target='%s.1'><data key='invisible'>1</data></edge>\n" % motif)
				for i in range (1,minVersion):
					f.write("<node id='%s.%i'><data key='invisible'>1</data></node>\n" %
							 (motif, i) )
					f.write("<edge source='%s.%i' target='%s.%i'><data key='invisible'>1</data></edge>\n" %
							 (motif,i,motif,i+1) )
			else:
				f.write("<edge source='root' target='%s.%i'><data key='invisible'>1</data></edge>\n" %
						(motif,minVersion) )		
		
			for i in range(minVersion,maxVersion):
				f.write("<edge source='%s.%i' target='%s.%i'></edge>\n" %
						(motif, i, motif, i+1) )
		else:
			f.write("<edge source='root' target='%s.%i'><data key='invisible'>1</data></edge>\n" %
					(motif, maxVersion) )
	
	for motif, parents in P.iteritems():
		t = motif.split('.')
		v = int(t[1])
		for parent in parents:
			f.write("<edge source='%s.%i' target='%s'></edge>\n" %
					 (parent, v-1, motif) )
			
	f.write('</graph>\n</graphml>')		
	f.close()
	
	# 	pdb.set_trace()
		
	db.close()
	print 'Done. Created file %s' % FILE


if __name__ == '__main__':
    aMotifsToGraphML()
    
        