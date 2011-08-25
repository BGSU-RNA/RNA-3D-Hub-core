"""

8/25/2011 Anton Petrov
Collection of common functions used in the Motif uploader

"""


def find_indices(tlist,m): 
	"""
	"""
	return [i for i,x in enumerate(tlist) if x == m]				


def transform_data(loops, groups):
	"""
	"""	
	return dict([(i, [loops[x] for x in find_indices(groups,i)]) \
				for i in set(groups)])
	





