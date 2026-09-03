
import pytest
import numpy as np

from qlsc import QLSC, GeoQLSC
from qlsc import q3c

# The integer grid coordinates (xi,yi) are the quantity the scheme is actually
# built from: an ipix value is the bit interleaving of a face number and a pair of
# integers in [0,nside). QLSC.ipix2xiyi decodes that (with the vendored xbits1 /
# ybits1 tables) and QLSC.xiyi2ipix encodes it (with the vendored q3c_xiyi2ipix and
# the xbits / ybits tables). The two are exact inverses of each other, which is
# what these tests pin: no trigonometry, no floating point, no tolerance anywhere.
#
# This matters because the same trip through (ra,dec) - xy2ang then ang2ipix - is
# only approximately invertible, and is at its worst exactly at a pixel boundary,
# which is where anything walking the grid spends all its time.

# Depths small enough to enumerate every pixel of the scheme (6 * 4**depth).
EXHAUSTIVE_DEPTHS = [0, 1, 2, 3, 4, 5]

# Depths too large to enumerate, sampled below. 30 is the maximum the scheme
# supports and the one where a float round trip loses; 26 is where an ipix value
# stops fitting exactly in a double.
SAMPLED_DEPTHS = [6, 10, 15, 20, 25, 26, 27, 29, 30]

ALL_DEPTHS = EXHAUSTIVE_DEPTHS + SAMPLED_DEPTHS


def _corner_and_edge_grid_positions(nside):
	'''
	The grid positions worth testing on one face: the four corners, the middle of
	each edge, one pixel in from each corner, and the middle of the face.
	'''
	last = nside - 1
	middle = nside // 2
	inward = min(1, last)
	values = sorted({0, inward, middle, last - inward, last})
	return [(xi, yi) for xi in values for yi in values]


@pytest.mark.parametrize("depth", EXHAUSTIVE_DEPTHS)
def test_round_trip_is_exact_for_every_pixel(depth):
	'''
	Test ipix -> (facenum,xi,yi) -> ipix over every pixel of the scheme, for every
	depth whose scheme is small enough to enumerate.
	'''
	q = QLSC(depth=depth)

	for ipix in range(q.nbins):
		facenum, xi, yi = q.ipix2xiyi(ipix)
		assert 0 <= facenum <= 5
		assert 0 <= xi < q.nside
		assert 0 <= yi < q.nside
		assert q.xiyi2ipix(facenum, xi, yi) == ipix


@pytest.mark.parametrize("depth", SAMPLED_DEPTHS)
def test_round_trip_is_exact_at_the_face_boundaries_and_corners(depth):
	'''
	Test the same round trip at the positions where a float round trip through the
	face plane is least trustworthy: the four corners of every face, the middles of
	its edges, and the pixels one step inside them.

	The first and last pixel of every face are included explicitly. At depth 30 the
	last pixel of the scheme is 6 * (2**30)**2 - 1 = 6917529027641081855, which is
	well past the 53 bits a double holds exactly - the same trap that once reported
	it as being on face 6.
	'''
	q = QLSC(depth=depth)

	for facenum in range(6):
		for xi, yi in _corner_and_edge_grid_positions(q.nside):
			ipix = q.xiyi2ipix(facenum, xi, yi)
			assert q.ipix2face(ipix) == facenum
			assert q.ipix2xiyi(ipix) == (facenum, xi, yi)

	first = 0
	last = q.nbins - 1
	assert q.ipix2xiyi(first) == (0, 0, 0)
	assert q.xiyi2ipix(*q.ipix2xiyi(last)) == last
	assert q.ipix2xiyi(last)[0] == 5


@pytest.mark.parametrize("depth", ALL_DEPTHS)
def test_round_trip_is_exact_for_sampled_pixels(depth):
	'''
	Test the round trip over a deterministic pseudo-random sample of the scheme, so
	that the interior of a face is covered too and not only its edges.
	'''
	q = QLSC(depth=depth)
	rng = np.random.default_rng(20260903 + depth)
	sample = rng.integers(0, q.nbins, size=min(500, q.nbins))

	for ipix in sample.tolist():
		assert q.xiyi2ipix(*q.ipix2xiyi(ipix)) == ipix


@pytest.mark.parametrize("depth", [0, 1, 4, 8, 16, 30])
def test_grid_coordinates_agree_with_ipix2xy(depth):
	'''
	Test that the integer grid position and the face-plane position of the same
	pixel are the same position: x is exactly (xi / nside) * 2 - 1, the rescaling
	the C decoder applies to the very same integers.

	The comparison is exact, not to a tolerance: both come from the same decode.
	'''
	q = QLSC(depth=depth)
	rng = np.random.default_rng(4242 + depth)
	sample = rng.integers(0, q.nbins, size=min(300, q.nbins))

	for ipix in sample.tolist():
		facenum_xy, x, y = q.ipix2xy(ipix)
		facenum_grid, xi, yi = q.ipix2xiyi(ipix)

		assert facenum_grid == facenum_xy
		assert x == (xi / q.nside) * 2 - 1
		assert y == (yi / q.nside) * 2 - 1


@pytest.mark.parametrize("depth", [2, 6, 12, 30])
def test_encode_agrees_with_ang2ipix(depth):
	'''
	Test that encoding the grid position of the pixel containing a point gives back
	the pixel ang2ipix names for that point. This ties the integer path to the
	scheme's own definition of where a point lands, rather than only to its inverse.
	'''
	q = QLSC(depth=depth)
	points = [(0., 0.), (10., 20.), (45., 35.26438968), (180., -45.), (270., 89.9),
			  (359.9, -89.9), (135., 0.), (315., 44.9)]

	for ra, dec in points:
		ipix = q.ang2ipix(ra, dec)
		assert q.xiyi2ipix(*q.ipix2xiyi(ipix)) == ipix


@pytest.mark.parametrize("depth", [2, 5, 30])
def test_stepping_the_grid_stays_inside_the_face(depth):
	'''
	Test that stepping xi and yi by one lands on a different pixel of the same face,
	and that the step is exact: the pixel it lands on decodes back to the stepped
	grid position.
	'''
	# nside >= 4 for every depth listed, so the middle pixel has a neighbour on
	# every side without leaving the face; the face-crossing cases belong to
	# ipix_neighbors, not here
	q = QLSC(depth=depth)
	facenum, xi, yi = 3, q.nside // 2, q.nside // 2
	here = q.xiyi2ipix(facenum, xi, yi)

	for dxi, dyi in [(1,0), (-1,0), (0,1), (0,-1), (1,1), (-1,-1), (1,-1), (-1,1)]:
		there = q.xiyi2ipix(facenum, xi + dxi, yi + dyi)
		assert there != here
		assert q.ipix2xiyi(there) == (facenum, xi + dxi, yi + dyi)
		assert q.ipix2face(there) == facenum


# ------------------------------------------------------------------ bad input

def test_xiyi2ipix_refuses_positions_off_the_face():
	'''
	Test that a grid coordinate outside [0,nside) is refused.

	It is not a wrapped position on this face: it belongs to a neighbouring face,
	whose grid is rotated with respect to this one. Reducing it modulo nside - which
	is what the vendored q3c_xiyi2ipix does with its table lookup - would answer with
	a pixel on the wrong side of the face, and a negative value indexes the
	interleaving table before its start.
	'''
	q = QLSC(depth=5)

	for xi, yi in [(-1, 0), (0, -1), (q.nside, 0), (0, q.nside), (-1, -1),
				   (q.nside + 1, 4), (4, q.nside + 1)]:
		with pytest.raises(ValueError):
			q.xiyi2ipix(2, xi, yi)


def test_xiyi2ipix_refuses_a_face_that_is_not_a_face():
	q = QLSC(depth=5)

	for facenum in [-1, 6, 7, 100]:
		with pytest.raises(ValueError):
			q.xiyi2ipix(facenum, 1, 1)


def test_xiyi2ipix_refuses_non_integers():
	'''
	Test that a fractional or boolean grid coordinate is refused rather than
	truncated. bool is a subclass of int in Python, and True is a mistake here, not
	the pixel 1.
	'''
	q = QLSC(depth=5)

	for bad in [1.5, 1.0, "1", None, True, False, np.float64(1.0)]:
		with pytest.raises(ValueError):
			q.xiyi2ipix(2, bad, 1)
		with pytest.raises(ValueError):
			q.xiyi2ipix(2, 1, bad)

	with pytest.raises(ValueError):
		q.xiyi2ipix(True, 1, 1)


def test_xiyi2ipix_requires_all_three_arguments():
	q = QLSC(depth=5)

	with pytest.raises(ValueError):
		q.xiyi2ipix()
	with pytest.raises(ValueError):
		q.xiyi2ipix(facenum=1)
	with pytest.raises(ValueError):
		q.xiyi2ipix(facenum=1, xi=1)


def test_xiyi2ipix_accepts_numpy_integers():
	'''
	Test that NumPy integers are accepted: they are what indexing or arithmetic on an
	ipix array produces, so refusing them would make the method awkward exactly where
	it is useful.
	'''
	q = QLSC(depth=5)
	expected = q.xiyi2ipix(2, 3, 4)

	assert q.xiyi2ipix(np.int64(2), np.int32(3), np.int8(4)) == expected


def test_ipix2xiyi_refuses_an_out_of_range_pixel():
	q = QLSC(depth=4)

	with pytest.raises(ValueError):
		q.ipix2xiyi(-1)
	with pytest.raises(ValueError):
		q.ipix2xiyi(q.nbins)


def test_ipix2xiyi_refuses_an_array():
	'''
	Test that array input is refused rather than answered for element zero. Unlike
	ipix2xy this method is not vectorized, and a silent element-zero read would be a
	wrong answer wearing the shape of a right one.
	'''
	q = QLSC(depth=4)

	for bad in [[1, 2], (1, 2), np.array([1, 2])]:
		with pytest.raises((TypeError, ValueError)):
			q.ipix2xiyi(bad)


def test_ipix2xiyi_refuses_a_fractional_pixel():
	q = QLSC(depth=4)

	with pytest.raises((TypeError, ValueError)):
		q.ipix2xiyi(1.5)


# ------------------------------------------------------------ the C functions

def test_c_level_xiyi2ipix_validates_its_arguments():
	'''
	Test that the C function refuses out-of-range values on its own, not only
	through the Python method: it is a public entry point of the extension module.
	'''
	q = QLSC(depth=5)

	with pytest.raises(ValueError):
		q3c.xiyi2ipix(q._hprm, -1, 0, 0)
	with pytest.raises(ValueError):
		q3c.xiyi2ipix(q._hprm, 6, 0, 0)
	with pytest.raises(ValueError):
		q3c.xiyi2ipix(q._hprm, 0, -1, 0)
	with pytest.raises(ValueError):
		q3c.xiyi2ipix(q._hprm, 0, 0, q.nside)
	with pytest.raises(TypeError):
		q3c.xiyi2ipix(q._hprm, 0, 1.5, 0)


def test_c_level_ipix2xiyi_validates_its_argument():
	q = QLSC(depth=5)

	with pytest.raises(ValueError):
		q3c.ipix2xiyi(q._hprm, q.nbins)
	with pytest.raises(TypeError):
		q3c.ipix2xiyi(q._hprm, np.array([1, 2]))


# ------------------------------------------------------------------ geographic

def test_the_geographic_class_inherits_the_grid_methods():
	'''
	Test that GeoQLSC gets both methods unchanged. Neither takes or returns a
	coordinate, so there is no longitude convention for the subclass to apply.
	'''
	g = GeoQLSC(depth=6)
	ipix = g.lonlat2ipix(lon=-74.0060, lat=40.7128)

	assert g.xiyi2ipix(*g.ipix2xiyi(ipix)) == ipix


# --------------------------------------------- the method the grid path fixed

@pytest.mark.parametrize("depth", [0, 1, 2, 5, 10, 20, 26, 30])
def test_xy2ipix_returns_the_pixel_ipix2xy_came_from(depth):
	'''
	Test that xy2ipix and ipix2xy are inverses at every depth.

	They were not. xy2ipix used to project the position back onto the sphere with
	xy2ang and hand the result to ang2ipix - two trips through trigonometry - and
	the corner ipix2xy returns sits exactly on a pixel boundary, where the last few
	bits decide the answer. At depth 30 that returned a neighbouring pixel for
	roughly 40% of the pixels of the scheme (1206 of a 3000-pixel sample). The
	integer grid path has no such boundary case: nside is a power of two, so a
	pixel's own corner scales back to its own grid coordinates exactly.
	'''
	q = QLSC(depth=depth)
	rng = np.random.default_rng(777 + depth)
	sample = (range(q.nbins) if q.nbins <= 500
			  else rng.integers(0, q.nbins, size=500).tolist())

	for ipix in sample:
		facenum, x, y = q.ipix2xy(ipix)
		assert q.xy2ipix(facenum=facenum, x=x, y=y) == ipix


@pytest.mark.parametrize("depth", [5, 20])
def test_xy2ipix_agrees_with_the_route_through_the_sky(depth):
	'''
	Test that positions inside a pixel - where the old route was never in doubt -
	still get the same answer as projecting onto the sphere and calling ang2ipix.

	This is the property the change had to preserve: the grid path is a different
	arithmetic for the same question, not a different question.
	'''
	q = QLSC(depth=depth)
	rng = np.random.default_rng(31337 + depth)

	for _ in range(400):
		facenum = int(rng.integers(0, 6))
		x = float(rng.uniform(-0.999, 0.999))
		y = float(rng.uniform(-0.999, 0.999))

		ra, dec = q.xy2ang(facenum=facenum, x=x, y=y)
		assert q.xy2ipix(facenum=facenum, x=x, y=y) == q.ang2ipix(ra, dec)


@pytest.mark.parametrize("depth", [1, 4, 30])
def test_xy2ipix_at_the_far_edge_of_a_face(depth):
	'''
	Test that x or y exactly at +1 - the far edge of the face, which scales to grid
	coordinate nside, one past the last pixel - lands on the last pixel rather than
	off the face. This is the adjustment q3c_ang2ipix_xy makes for the same reason.
	'''
	q = QLSC(depth=depth)
	last = q.nside - 1

	assert q.xy2ipix(facenum=2, x=1., y=1.) == q.xiyi2ipix(2, last, last)
	assert q.xy2ipix(facenum=2, x=1., y=-1.) == q.xiyi2ipix(2, last, 0)
	assert q.xy2ipix(facenum=2, x=-1., y=1.) == q.xiyi2ipix(2, 0, last)
	assert q.xy2ipix(facenum=2, x=-1., y=-1.) == q.xiyi2ipix(2, 0, 0)


def test_xy2ipix_still_validates_its_arguments():
	'''
	Test that the validation xy2ipix used to inherit from xy2ang is still applied now
	that it no longer calls it.
	'''
	q = QLSC(depth=8)

	with pytest.raises(ValueError):
		q.xy2ipix(x=0.1, y=0.2)                  # no face
	with pytest.raises(ValueError):
		q.xy2ipix(facenum=6, x=0.1, y=0.2)       # not a face
	with pytest.raises(ValueError):
		q.xy2ipix(facenum=-1, x=0.1, y=0.2)
	with pytest.raises(ValueError):
		q.xy2ipix(facenum=0.5, x=0.1, y=0.2)     # a float where a face number belongs
	with pytest.raises(ValueError):
		q.xy2ipix(facenum=2, x=0.1)              # missing y
	with pytest.raises(ValueError):
		q.xy2ipix(facenum=2, y=0.2)              # missing x
	with pytest.raises(ValueError):
		q.xy2ipix(facenum=2, x=1.4, y=0.2)       # past the edge of the face
	with pytest.raises(ValueError):
		q.xy2ipix(facenum=2, x=0.1, y=-1.4)
	with pytest.raises((ValueError, TypeError)):
		q.xy2ipix(facenum=2, x=float("nan"), y=0.2)
