
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
		idx.radial_query(ra=12., dec=34., radius=0.1)

def test_radial_query_negative_radius():
	'''
	Test that a negative radius raises a ValueError.
	'''
	idx = QLSCIndex(qlsc=QLSC(depth=30))
	idx.add_point(12., 34.)
	with pytest.raises(ValueError):
		idx.radial_query(ra=12., dec=34., radius=-1.)

def test_radial_query_zero_radius():
	'''
	Test that a radius of zero returns an empty result rather than an SQL error.
	'''
	idx = QLSCIndex(qlsc=QLSC(depth=30))
	idx.add_point(12., 34.)
	matches = idx.radial_query(ra=12., dec=34., radius=0.)
	assert len(matches) == 0
	matches = idx.radial_query(ra=12., dec=34., radius=0., return_key=True)
	assert len(matches) == 0
