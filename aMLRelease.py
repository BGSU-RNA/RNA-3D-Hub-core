"""



"""


import csv
import functions as af

__author__ = 'Anton Petrov'

class Release:

	def __init__(self, Motifs = {}, loops = [], groups = [], history = {}, 
	             version = 0):
		self.M = Motifs  # dict( [motif]:[loop1..loopN] )
		self.l = loops
		self.g = groups
		self.history = history
		self.version = version
		self.message = ''
		self.newRelease = False


	def set_version(self, version):
		self.version = version


	def set_release(self, release):
		self.newRelease = release
		
		
	def set_message(self, message):
		self.message = message


	def load_from_file(self,file):
		"""
		"""
		self.l, self.g = self.read_motifs_csv_file(file)
		self.M = af.transform_data(self.l, self.g)
				

	def read_motifs_csv_file(self, file):
		"""
		"""
		motifReader = csv.reader(open(file), delimiter=',', quotechar='"')
		loops = []
		groups = []
		for row in motifReader:
			loops.append(row[0])
			groups.append(row[1])
		return loops, groups


	def load_from_database(self,db):
		"""
		"""
		self.l, self.g = db.get_current_release()
		self.M = af.transform_data(self.l, self.g)					
