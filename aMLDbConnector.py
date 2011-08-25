"""

8/24/2011 Anton Petrov
Class for interfacing with the MotifAtlas database

Public methods:
	.query('command')
	.generateId()
	.close()
	.load_motifs_data({motif_id:[loop1..loopN]}, version):
	.load_history_data({motif_id:[parent1..parentN]}, version)
	.get_current_version()
	.get_versions()	
	.get_current_release()
	
Public variables:
	.motifsTable
	.historyTable
	
"""

import sys
import pdb
import random
import MySQLdb as mdb

import aMLRelease as amr
import functions as af

# suppresses mysql warning triggered by IF NOT EXISTS query
from warnings import filterwarnings
import MySQLdb as Database
filterwarnings('ignore', category = Database.Warning)


# default initialization parameters
HOST	 = 'localhost'
USER	 = 'root'
PASSWORD = 'rumata7'
DATABASE = 'MotifAtlas1'
MOTIFS_TABLE   = 'Motifs'
HISTORY_TABLE  = 'History'
VERSIONS_TABLE = 'Versions'


class aDbConnector():

	def __init__(self, database = DATABASE, host = HOST, user = USER, 
	             password = PASSWORD, motifsTable = MOTIFS_TABLE, 
	             historyTable = HISTORY_TABLE, versionsTable = VERSIONS_TABLE):
		self.__database = database
		self.__host = host
		self.__user = user
		self.__password = password
		self.motifsTable   = motifsTable
		self.historyTable  = historyTable
		self.versionsTable = versionsTable
		# initialize database connection
 		self.__connect_to_mysql()
 		self.__select_database()
 		self.__create_all_tables()
 		
################################################################################
# Public functions
################################################################################

	def upload_first_release(self,R):
		""" R - instance of the Release class
		"""
		# some error checking
		if R.version == 0:
			print 'Tried to upload a release with version 0'
			sys.exit()
		if len(R.M) == 0:
			print 'Tried to upload a release with an empty motif structure'
			sys.exit()
		
		# generate ids for all loops 
		N = dict()
		for group, loops in R.M.iteritems():
			newId = self.__generate_id()
			N[newId] = loops
		
		R.M = N
		self.__insert_motifs_data(R)
		self.__insert_versions_data(R)

								
	def get_current_version(self):
		"""
		"""	
		self.execute("SELECT MAX(version) FROM %s;" % self.motifsTable);
		result = self.__fetchone()
		
		if result is None:
			Version = 0
		else:
			Version = result[0]
			
		print "Current version is %i"   % Version
		return Version
		
				
	def __generate_id(self):
		"""
		"""
		prefix = 'IL'
		while True:
			id = '%s%05d' % (prefix, random.randint(1, 99999) )
			self.execute("SELECT motif_id FROM %s WHERE motif_id='%s';" 
			             % (self.motifsTable, id) )
			result = self.__fetchone()
			if result is None: break
		
		return id
		
		
################################################################################		
# Private functions						
# Insert data into tables
################################################################################

	def __insert_motifs_data(self, R):
		""" loop_id, motif_id, version_id
		"""		
		for motif, loops in R.M.iteritems():
			for loop in loops:
				command = "INSERT INTO %s VALUES ('%s','%s',%s);" \
				           % (self.motifsTable, loop, motif, R.version) 
				self.execute(command)

		self.conn.commit()	
		print 'Motifs data loaded'
		
	
	def __insert_history_data(self, M, version):
		""" motif_id, parents, version_id
		"""		
		for child, parents in M.iteritems():
			concat = reduce(lambda x, y: x+','+y, parents)
			concat = concat[:-1]

			command = "INSERT INTO %s VALUES ('%s','%s',%s);" \
					   % (self.historyTable, child,concat,version)
			self.execute(command)
		
		self.conn.commit()	
		print 'History data loaded'

						
	def __insert_versions_data(self, R):
		"""
		"""
		# get the current release_id
		command = "SELECT MAX(release_id) FROM %s;" % self.versionsTable
		self.execute(command)
		result = self.__fetchone()
		
		# set to 1 for the first time, increment if R.newRelease is 1
		if result[0] is None:
			release = 1
		elif R.newRelease == True:
			release = result[0] + 1		
		else:
			release = result[0]

		pdb.set_trace()				
		command = "INSERT INTO %s (release_id, version, comments) \
		           VALUES (%i, %i, '%s');" % ( self.versionsTable, release, R.version, R.message )
		self.execute(command)
        
        print 'Updated versions table'

	
################################################################################	
# Manipulate existing releases
################################################################################

	def get_all_releases(self):
		pass						
	
	
	def get_last_release(self):
		"""
		"""
		self.execute("SELECT MAX(version) FROM %s;" % self.motifsTable) 
		version = self.__fetchone()
		
		if version[0] is not None:
			return self.get_release(version[0])
		else:
			return amr.Release()
			
							
	def get_release(self, version):
		"""
		"""
		self.cursor.execute("SELECT loop_id, motif_id FROM %s WHERE version = %i;"\
							% (self.motifsTable, version) )
		rows = self.cursor.fetchall()
		l = []; g = [];
		for row in rows:
			l.append(row[0])
			g.append(row[1])
		
		M = af.transform_data(l, g)
		
		R = amr.Release(loops = l, groups = g, Motifs = M, version = version)
		
		return R
		
	
	def delete_release(self, version):
		"""
		"""
		command = "DELETE * FROM %s WHERE version = %i;" % self.motifsTable
		self.execute(command)
	
	
################################################################################
# Mysql connection and table management
################################################################################	

	def __connect_to_mysql(self):
		""" Establish database connection.
		"""
		try:
			self.conn = mdb.connect(self.__host, self.__user, self.__password)
			self.cursor = self.conn.cursor()
		except mdb.Error, e:  
			print "Error %d: %s" % (e.args[0],e.args[1])
 			sys.exit(1)

		
	def __select_database(self):
		""" Create and use the database.
		"""
		try:	
			self.cursor.execute("CREATE DATABASE IF NOT EXISTS %s;" % self.__database)
			self.cursor.execute("USE %s;" % self.__database)	
		except mdb.Error, e:
			print "Error %d: %s" % (e.args[0],e.args[1])
			sys.exit(1)	
	
	
	def __create_all_tables(self):
		"""
		"""
		self.__create_motifs_table()
		self.__create_history_table()		
		self.__create_versions_table()
	
	
	def __create_motifs_table(self):
		""" loop_id, motif_id, version
		"""
		try:	
			self.cursor.execute("""CREATE TABLE IF NOT EXISTS %s \
			                    (loop_id VARCHAR(11), motif_id VARCHAR(11),\
			                    version INT);""" % self.motifsTable)
		except mdb.Error,e:
			print "Error %d: %s" % (e.args[0],e.args[1])
	
		
	def __create_history_table(self):
		""" motif_id, parents, version
		"""
		try:	
			self.cursor.execute("""CREATE TABLE IF NOT EXISTS %s \
			                    (motif_id VARCHAR(11), parents TEXT,\
			                    version INT);""" % self.historyTable)
		except mdb.Error,e:
			print "Error %d: %s" % (e.args[0],e.args[1])
	

	def __create_versions_table(self):
		""" id, release, version, comments
		"""
		try:	
			self.cursor.execute("""CREATE TABLE IF NOT EXISTS %s
			                     (id MEDIUMINT NOT NULL AUTO_INCREMENT,
			                     release_id INT, version INT, comments TEXT,
			                     PRIMARY KEY (id) );""" 
			                     % self.versionsTable)
		except mdb.Error,e:
			print "Error %d: %s" % (e.args[0],e.args[1])		
	
	
################################################################################
# Housekeeping functions
################################################################################
					
	def execute(self,command):
		""" 
		""" 
		try:
			self.cursor.execute(command)
	
		except mdb.Error,e:
			print "Error %d: %s" % (e.args[0],e.args[1])
			sys.exit()


	def close(self):
		self.cursor.close()
		self.conn.close()


	def __fetchone(self):
		return self.cursor.fetchone()


	def __fetchall(self):
		return self.cursor.fetchall()

	
################################################################################


	
	# 	def get_versions(self):
# 		"""
# 		"""	
# 		self.cursor.execute("SELECT MAX(version) FROM %s;" % self.motifsTable);
# 		result = self.cursor.fetchone()
# 		
# 		if result[0] is None:
# 			newVersion = 1
# 			oldVersion = 0
# 		else:
# 			oldVersion = result[0]
# 			newVersion = oldVersion + 1
# 			
# 		print "Current version is %i"   % oldVersion
# 		print "Next version will be %s" % newVersion
# 		return newVersion, oldVersion

		
# 	def get_current_release(self):
# 		"""
# 		"""		
# 		self.cursor.execute("SELECT loop_id, motif_id FROM %s WHERE version = \
# 							(SELECT MAX(version) FROM %s);" 
# 						    % (self.motifsTable, self.motifsTable) )
# 		rows = self.cursor.fetchall()
# 		L = []; G = []
# 		for row in rows:
# 			L.append(row[0])
# 			G.append(row[1])
# 					
# 		return L, G				

