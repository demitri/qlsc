
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

# (10,20) and (0,0) are positions where q3c_radial_query produces no candidate
# ranges at all for tiny radii; (12,34) is a position where it does; the last
# position is one where a 1e-6 degree search radius overflows q3c's internal
# range capacity and the search must be retried with a larger radius
@pytest.mark.parametrize("ra, dec", [(10., 20.), (0., 0.), (12., 34.), (45., 45.), (200., -60.),
                                     (316.73891370708168, 3.0827793488343289)])
def test_radial_query_zero_radius(ra, dec):
	'''
	Test that a radius of zero returns only exactly coincident points
	(the radius comparison is inclusive), never an SQL error.
	'''
	idx = QLSCIndex(qlsc=QLSC(depth=30))
	idx.add_point(ra, dec)
	idx.add_point(ra + 1., dec)
	matches = idx.radial_query(ra=ra, dec=dec, radius=0.)
	assert matches.shape == (2,) # the single exact match, squeezed to (ra,dec)
	matches = idx.radial_query(ra=ra + 0.5, dec=dec, radius=0.)
	assert len(matches) == 0
	matches = idx.radial_query(ra=ra + 0.5, dec=dec, radius=0., return_key=True)
	assert len(matches) == 0

# stored representation vs equivalent query representation of the same point
@pytest.mark.parametrize("stored, query", [
	((360., 0.), (0., 0.)),     # ra seam
	((-360., 20.), (0., 20.)),  # ra seam, negative
	((180., 90.), (0., 90.)),   # north pole: ra is irrelevant
	((0., -90.), (37., -90.))   # south pole
])
def test_radial_query_zero_radius_equivalent_coordinates(stored, query):
	'''
	Test that a zero-radius query matches a stored point given in an
	equivalent coordinate representation (ra periodicity, poles).
	'''
	idx = QLSCIndex(qlsc=QLSC(depth=30))
	idx.add_point(*stored)
	matches = idx.radial_query(ra=query[0], dec=query[1], radius=0.)
	assert matches.shape == (2,), f"stored {stored} not found by zero-radius query at {query}"

@pytest.mark.parametrize("ra, dec", [(10., 20.), (0., 0.), (12., 34.)])
def test_radial_query_tiny_radius(ra, dec):
	'''
	Test tiny radii for which q3c_radial_query can produce no candidate ranges;
	the query must still find indexed points within the radius.
	'''
	idx = QLSCIndex(qlsc=QLSC(depth=30))
	idx.add_point(ra, dec)
	idx.add_point(ra + 1., dec)
	matches = idx.radial_query(ra=ra, dec=dec, radius=1e-8)
	assert matches.shape == (2,) # the coincident point only
