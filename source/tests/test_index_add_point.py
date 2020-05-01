
import pytest
import numpy as np

from qlsc import QLSC, QLSCIndex

def test_add_point_keyword_parameter():
	'''
	Test QLSC.add_point using keyword parameters.
	'''
	q30 = QLSC(depth=30)
	idx = QLSCIndex(qlsc=q30)
	
	ra = 12
	dec = 34
	
	idx.add_point(ra=ra, dec=dec)
	
	# make sure point is there
	matches = idx.radial_query(ra=ra, dec=dec, radius=10/3600)
	
	assert np.any(matches == [ra,dec])
	
def test_add_point_position_parameter():
	'''
	Test QLSC.add_point using position parameters.
	'''
	q30 = QLSC(depth=30)
	idx = QLSCIndex(qlsc=q30)
	
	ra = 12
	dec = 34
	
	idx.add_point(ra,dec)
	
	# make sure point is there
	matches = idx.radial_query(ra=ra, dec=dec, radius=10/3600)
	
	assert np.any(matches == [ra,dec])
	
def test_add_point_retrieve_key():
	'''
	Test that a key added with a point is returned.
	'''
	q30 = QLSC(depth=30)
	idx = QLSCIndex(qlsc=q30)
	
	ra = 12
	dec = 34
	key = "my key"
	
	idx.add_point(ra,dec,key)
	
	# make sure point is there
	matches = idx.radial_query(ra=ra, dec=dec, radius=10/3600, return_key=True)
	
	# find record, note that keys are not required to be unique
	record =  matches[matches['key'] == key][0] # returns a recarray, get the tuple with [0]
	
	assert record[2] == key, "a provided coordinate key was not correctly retrieved"
	
def test_add_point_duplicate_points_ignored():
	'''
	Test that duplicate ra,dec pairs are ignored.
	'''
	
	q30 = QLSC(depth=30)
	idx = QLSCIndex(qlsc=q30)
	
	ra = 12
	dec = 34
	key = "my key"
	
	idx.add_point(ra,dec,key)
	idx.add_point(ra,dec,key)
	idx.add_point(ra,dec,key)
	idx.add_point(ra,dec,key)
	idx.add_point(ra,dec,key)

	assert idx.number_of_points == 1, "duplicate points are being incorrectly added"
