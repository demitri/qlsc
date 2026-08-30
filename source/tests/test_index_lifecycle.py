
import pytest
import numpy as np

from qlsc import QLSC, QLSCIndex

def test_context_manager():
	'''
	Test that QLSCIndex works as a context manager.
	'''
	with QLSCIndex(qlsc=QLSC(depth=30)) as idx:
		idx.add_point(12., 34.)
		matches = idx.radial_query(ra=12., dec=34., radius=0.1)
		assert len(matches) == 2 # single match returned as squeezed (ra,dec)

def test_close_is_idempotent():
	'''
	Test that close() can be called more than once.
	'''
	idx = QLSCIndex(qlsc=QLSC(depth=30))
	idx.close()
	idx.close()

def test_use_after_close_raises():
	'''
	Test that using a closed index raises rather than silently
	operating on a fresh, empty database.
	'''
	idx = QLSCIndex(qlsc=QLSC(depth=30))
	idx.add_point(12., 34.)
	idx.close()
	with pytest.raises(RuntimeError):
		idx.add_point(56., 78.)
	with pytest.raises(RuntimeError):
		idx.add_points(ra=np.array([56.]), dec=np.array([78.]))
	with pytest.raises(RuntimeError):
		idx.radial_query(ra=12., dec=34., radius=0.1)
	with pytest.raises(RuntimeError):
		# a zero radius must not short-circuit the closed check
		idx.radial_query(ra=12., dec=34., radius=0.)
	with pytest.raises(RuntimeError):
		idx.number_of_points

def test_radial_query_radius_out_of_range():
	'''
	Test that radii outside [0,180] raise a ValueError.
	'''
	idx = QLSCIndex(qlsc=QLSC(depth=30))
	idx.add_point(12., 34.)
	with pytest.raises(ValueError):
		idx.radial_query(ra=12., dec=34., radius=-1.)
	with pytest.raises(ValueError):
		# above 180 degrees the sin^2(radius/2) exact filter becomes periodic
		idx.radial_query(ra=12., dec=34., radius=181.)

def test_radial_query_radius_inclusive():
	'''
	Test that the radius comparison is inclusive: a radius of 180 degrees
	returns the whole sphere, including an exactly antipodal point.
	'''
	idx = QLSCIndex(qlsc=QLSC(depth=30))
	idx.add_point(12., 34.)
	idx.add_point(192., -34.) # antipode of (12,34)
	matches = idx.radial_query(ra=12., dec=34., radius=180.)
	assert matches.shape == (2, 2)

def test_radial_query_zero_radius():
	'''
	Test that a radius of zero returns only exactly coincident points
	(the radius comparison is inclusive), never an SQL error.
	'''
	idx = QLSCIndex(qlsc=QLSC(depth=30))
	idx.add_point(12., 34.)
	idx.add_point(13., 34.)
	matches = idx.radial_query(ra=12., dec=34., radius=0.)
	assert matches.shape == (2,) # the single exact match, squeezed to (ra,dec)
	matches = idx.radial_query(ra=50., dec=50., radius=0.)
	assert len(matches) == 0
	matches = idx.radial_query(ra=50., dec=50., radius=0., return_key=True)
	assert len(matches) == 0
