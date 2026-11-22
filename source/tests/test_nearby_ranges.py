
import pytest
import numpy as np

from qlsc import QLSC, GeoQLSC
from qlsc import q3c
from qlsc.generate import sunflower_points_on_sphere

# QLSC.nearby_ipix_ranges exposes q3c_get_nearby: Q3C's cheap "which squares of the
# quadtree are near this cone" calculation, the one its PostgreSQL extension uses to
# seed a positional join. It does no refinement and emits at most four ranges, each
# a single quadtree square at the level just coarse enough to hold the cone's
# bounding box on a cube face.
#
# READ THIS BEFORE CHANGING THE COVERAGE TEST BELOW. The method's docstring says
# the ranges are not a guaranteed cover, and that is deliberate: on a cube face the
# cone spills onto, q3c emits exactly one square. Before q3c v2.0.6 that square sat
# at a corner of the face's bounding box with nothing bounding the rest of the box
# to it (upstream found ellipses it missed; fixed in v2.0.6, which grows the square
# up the quadtree until it holds the whole box). test_nearby_ranges_cover_the_cone
# pins coverage over a deterministic grid. It is a measurement, not a proof of an
# invariant - if it ever fails, the finding is real and belongs upstream, not in a
# widened tolerance here.


def _angular_distance(ra0, dec0, ra, dec):
	'''
	Angular separation in degrees between one point and an array of them, computed
	independently of qlsc's own distance functions so that the coverage test does
	not check q3c against itself.
	'''
	r0, d0 = np.radians(ra0), np.radians(dec0)
	r, d = np.radians(ra), np.radians(dec)
	cosine = np.sin(d0) * np.sin(d) + np.cos(d0) * np.cos(d) * np.cos(r - r0)
	return np.degrees(np.arccos(np.clip(cosine, -1., 1.)))


def _covered(ranges, ipix_values):
	'''
	A boolean array saying which of the given ipix values fall in one of the
	inclusive [start,end] ranges.
	'''
	covered = np.zeros(np.shape(ipix_values), dtype=bool)
	for start, end in ranges:
		covered |= (ipix_values >= start) & (ipix_values <= end)
	return covered


# Centres chosen to sit on the features of the cube: face centres, face boundaries
# (ra = 45 + 90k), the poles, and the latitude of the cube vertices (35.26438968 deg,
# atan(1/sqrt(2)) in degrees), where three faces meet.
CENTRES = [(0., 0.), (45., 0.), (10., 20.), (90., 45.), (135., -44.),
		   (180., 35.26438968), (225., -35.26438968), (270., 89.9),
		   (315., -89.9), (359.9, 1.), (45., 35.26438968), (30., -20.)]

RADII = [0., 1e-7, 1e-3, 0.1, 1., 5., 20., 34.9]


@pytest.mark.parametrize("depth", [2, 5, 8])
@pytest.mark.parametrize("centre", CENTRES)
def test_nearby_ranges_cover_the_cone(depth, centre):
	'''
	Test that every pixel holding a point inside the cone is inside the returned
	ranges, over a deterministic grid of centres, radii and depths.

	This is an empirical pin, not a proof - see the note at the top of this file and
	the warning in the method's docstring.
	'''
	q = QLSC(depth=depth)
	points = sunflower_points_on_sphere(n=50000)
	ipix = q.ang2ipix(points[:,0], points[:,1])
	ra0, dec0 = centre

	for radius in RADII:
		separation = _angular_distance(ra0, dec0, points[:,0], points[:,1])
		inside = separation <= radius
		if not inside.any():
			continue

		ranges = q.nearby_ipix_ranges(ra0, dec0, radius)
		missed = int((~_covered(ranges, ipix[inside])).sum())

		assert missed == 0, (f"{missed} of {int(inside.sum())} points inside a cone at "
							 f"(ra={ra0}, dec={dec0}) with radius={radius} deg fell outside the "
							 f"{len(ranges)} nearby range(s) at depth {depth}")


@pytest.mark.parametrize("depth", [0, 1, 5, 20, 30])
def test_the_ranges_are_well_formed(depth):
	'''
	Test the shape and content of the result: an (n,2) int64 array of at most four
	ranges, each one running forwards, every value a pixel of this scheme.
	'''
	q = QLSC(depth=depth)

	for ra, dec in CENTRES:
		for radius in RADII + [35., 90., 180.]:
			ranges = q.nearby_ipix_ranges(ra, dec, radius)

			assert ranges.dtype == np.int64
			assert ranges.ndim == 2 and ranges.shape[1] == 2
			assert 1 <= ranges.shape[0] <= 4

			assert (ranges[:,0] <= ranges[:,1]).all()
			assert (ranges >= 0).all()
			assert (ranges < q.nbins).all()


@pytest.mark.parametrize("depth", [1, 6, 30])
def test_a_radius_past_q3c_maxrad_is_the_whole_sky(depth):
	'''
	Test that past Q3C_MAXRAD (35 degrees) q3c stops doing geometry and answers with
	the whole sky, and that the sentinel it spells [-1, 6*nside**2] is reported here
	as the pixels of this scheme.
	'''
	q = QLSC(depth=depth)

	for radius in [35.001, 45., 90., 180., 360.]:
		ranges = q.nearby_ipix_ranges(10., 20., radius)
		assert ranges.tolist() == [[0, q.nbins - 1]]

	# the C function's own spelling of the same answer, for the record
	raw = q3c.nearby_query(q._hprm, 10., 20., 90.)
	assert raw.tolist() == [[-1, 6 * q.nside * q.nside]]


def test_a_zero_radius_names_the_pixel_containing_the_centre():
	'''
	Test that a cone of no size still returns a range, and that it contains the pixel
	the centre itself falls in.
	'''
	q = QLSC(depth=10)

	for ra, dec in CENTRES:
		ranges = q.nearby_ipix_ranges(ra, dec, 0.)
		assert _covered(ranges, np.array([q.ang2ipix(ra, dec)])).all()


@pytest.mark.parametrize("centre", CENTRES)
def test_the_ranges_are_periodic_in_right_ascension(centre):
	'''
	Test that the same cone written with a shifted right ascension gives the same
	ranges. q3c_get_nearby reaches q3c_get_facenum (through q3c_get_region_facenum),
	which is documented "strictly 0<=ra<360" and does not enforce it, so the wrapper
	reduces ra before the call - the same class of defect this repository's
	test_ra_wrapping.py was written for.
	'''
	q = QLSC(depth=12)
	ra, dec = centre

	for radius in [0.01, 1., 20.]:
		expected = q.nearby_ipix_ranges(ra, dec, radius).tolist()
		for shift in [-1080, -720, -360, 360, 720, 1080]:
			assert q.nearby_ipix_ranges(ra + shift, dec, radius).tolist() == expected


def test_the_ranges_are_coarser_than_a_radial_query():
	'''
	Test the property that is the whole point of the method: it is cheap and coarse.
	The refined ranges of a cone search cover fewer pixels than the nearby squares
	that contain them.

	This is what makes it a pre-filter and not an answer.
	'''
	q = QLSC(depth=20)
	ra, dec, radius = 10., 20., 0.5

	nearby_size = int(sum(end - start + 1 for start, end in q.nearby_ipix_ranges(ra, dec, radius)))

	fulls, partials = q3c.radial_query(q._hprm, ra, dec, radius)
	query_size = int(sum(end - start + 1 for start, end in np.vstack([fulls, partials])))

	assert query_size < nearby_size


# ------------------------------------------------------------------ bad input

def test_a_negative_radius_is_refused():
	'''
	Test that a negative radius raises rather than returning an empty answer.

	q3c_get_nearby answers a negative radius with four filler pairs, which arrive
	here looking exactly like "nothing is nearby" - a wrong answer that a caller has
	no way to tell from a right one.
	'''
	q = QLSC(depth=8)

	with pytest.raises(ValueError):
		q.nearby_ipix_ranges(10., 20., -1.)
	with pytest.raises(ValueError):
		q3c.nearby_query(q._hprm, 10., 20., -1.)


def test_non_finite_values_are_refused():
	q = QLSC(depth=8)

	for bad in [float("nan"), float("inf"), float("-inf")]:
		with pytest.raises(ValueError):
			q.nearby_ipix_ranges(bad, 20., 1.)
		with pytest.raises(ValueError):
			q.nearby_ipix_ranges(10., bad, 1.)
		with pytest.raises(ValueError):
			q.nearby_ipix_ranges(10., 20., bad)
		with pytest.raises(ValueError):
			q3c.nearby_query(q._hprm, 10., bad, 1.)


def test_a_declination_out_of_range_is_refused_at_the_c_level():
	'''
	Test the C function's own declination check. The Python method normalizes an
	out-of-range declination instead (a pole crossing is a real position), which is
	what the next test covers.
	'''
	q = QLSC(depth=8)

	with pytest.raises(ValueError):
		q3c.nearby_query(q._hprm, 10., 91., 1.)


def test_a_declination_past_the_pole_is_normalized():
	'''
	Test that a declination past the pole is rewritten as the position it names,
	matching every other coordinate-taking method in the package.
	'''
	q = QLSC(depth=8)

	assert (q.nearby_ipix_ranges(10., 95., 1.).tolist() ==
			q.nearby_ipix_ranges(190., 85., 1.).tolist())


def test_all_three_arguments_are_required():
	q = QLSC(depth=8)

	with pytest.raises(ValueError):
		q.nearby_ipix_ranges()
	with pytest.raises(ValueError):
		q.nearby_ipix_ranges(10., 20.)
	with pytest.raises(ValueError):
		q.nearby_ipix_ranges(ra=10., radius=1.)


# ------------------------------------------------------------------ geographic

def test_the_geographic_spelling():
	'''
	Test that GeoQLSC accepts lon,lat and answers identically. Pixel numbers carry no
	longitude convention, so the result is the inherited one unchanged.
	'''
	g = GeoQLSC(depth=12)

	assert (g.nearby_ipix_ranges(lon=-74.0060, lat=40.7128, radius=1.).tolist() ==
			g.nearby_ipix_ranges(ra=285.9940, dec=40.7128, radius=1.).tolist())
