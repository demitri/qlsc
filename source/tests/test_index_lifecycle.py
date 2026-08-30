
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

def test_number_of_points_empty_index():
	'''
	Test that number_of_points is 0 for an empty index (not None).
	'''
	idx = QLSCIndex(qlsc=QLSC(depth=30))
	assert idx.number_of_points == 0
	idx.add_point(12., 34.)
	assert idx.number_of_points == 1

def test_add_points_single_point():
	'''
	Test add_points with a single point in all three argument forms.

	Regression test: ang2ipix collapsed a length-1 array input to a bare int
	(a size-1 array satisfies the C extension's scalar signature), so
	add_points raised TypeError when given exactly one point.
	'''
	for kwargs in ({"points": np.array([[10., 30.]])},
	               {"ra": np.array([10.]), "dec": np.array([30.])},
	               {"points": np.array([[10., 100.]])}): # out-of-range dec, n=1
		idx = QLSCIndex(qlsc=QLSC(depth=30))
		idx.add_points(**kwargs)
		assert idx.number_of_points == 1

def test_ang2ipix_length_1_array():
	'''
	Test that ang2ipix returns an array for array input, even length-1
	(per its documented contract).
	'''
	q = QLSC(depth=30)
	result = q.ang2ipix(ra=np.array([10.]), dec=np.array([30.]))
	assert isinstance(result, np.ndarray)
	assert result.shape == (1,)
	assert result[0] == q.ang2ipix(10., 30.)

def test_add_points_normalizes_stored_coordinates():
	'''
	Test that add_points normalizes out-of-range coordinates both for the ipix
	calculation and for the STORED ra,dec values.

	Regression test: the normalized values were previously computed and discarded,
	storing e.g. dec=100 verbatim; and the ra=/dec= argument form raised a spurious
	ValueError for out-of-range input.
	'''
	for kwargs in ({"points": np.array([[10., 100.], [20., 30.]])},
	               {"ra": np.array([10., 20.]), "dec": np.array([100., 30.])}):
		idx = QLSCIndex(qlsc=QLSC(depth=30))
		idx.add_points(**kwargs)
		# (10,100) normalizes to (190,80)
		matches = idx.radial_query(ra=190., dec=80., radius=0.5)
		assert matches.shape == (2,)
		ra, dec = matches
		assert -90 <= dec <= 90, f"stored dec value {dec} was not normalized"
		assert (ra, dec) == (190., 80.)

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
