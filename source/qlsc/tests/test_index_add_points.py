
import pytest
import numpy as np

from qlsc import QLSC, QLSCIndex
from qlsc.generate import sunflower_points_on_sphere

def test_add_points_radec_keyword_parameter():
	'''
	Test QLSC.add_points using ra,dec and keyword parameters.
	'''
	q30 = QLSC(depth=30)
	idx = QLSCIndex(qlsc=q30)
	
	sample_points = sunflower_points_on_sphere(n=100)
	
	ra = sample_points[:,0]
	dec = sample_points[:,1]
	
	idx.add_points(ra=ra, dec=dec)
		
	assert idx.number_of_points == 100, "points don't appear to be in the database"
	
def test_add_points_with_single_array():
	'''
	Test QLSC.add_point using position parameters.
	'''
	q30 = QLSC(depth=30)
	idx = QLSCIndex(qlsc=q30)
	
	sample_points = sunflower_points_on_sphere(n=100)
	
	idx.add_points(points=sample_points)
	
	assert idx.number_of_points == 100, "points don't appear to be in the database"
		
def test_add_points_duplicate_points_ignored():
	'''
	Test that duplicate ra,dec pairs are ignored.
	'''
	
	q30 = QLSC(depth=30)
	idx = QLSCIndex(qlsc=q30)
	
	sample_points = sunflower_points_on_sphere(n=100)
	
	idx.add_points(points=sample_points)
	idx.add_points(points=sample_points)
	idx.add_points(points=sample_points)
	idx.add_points(points=sample_points)
	idx.add_points(points=sample_points)

	assert idx.number_of_points == 100, "duplicate points are being incorrectly added"

# could add testing with keys...
