
import pytest
import numpy as np

from qlsc import QLSC, q3c
from qlsc.qlsc_base import ELLIPSE_MAXIMUM_SEMI_MAJOR_AXIS

# Q3C stops doing geometry for a region of Q3C_MAXRAD degrees or more and returns the
# whole sphere as one partially covered range: above that size a cone or an ellipse can
# project to a parabola or a hyperbola on a cube face and the quadtree machinery breaks.
#
# That fallback is load-bearing for qlsc, not an implementation detail of q3cube.c. Every
# region query here rests on one invariant - the coarse ranges are a SUPERSET of the
# region - and for any query at or above Q3C_MAXRAD the fallback is the only thing
# providing it. If an upstream sync ever removed or narrowed the fallback, large queries
# would start under-covering and silently lose points, with no test failing anywhere else:
# the sizes involved are larger than the ones the query tests use.
#
# These tests pin the fallback itself, at the C boundary, so that cannot happen quietly.
# scripts/sync_q3c.py runs the full suite after a sync.

q30 = QLSC(depth=30)

def whole_sphere_partials(qlsc:QLSC) -> list:
	'''
	The single ipix range q3c emits for a whole-sphere query, at the given scheme.

	The lower bound is -1 rather than 0: q3c writes a half-open range that starts below
	the first pixel, and the Python layer's range mapping copes with it.
	'''
	return [[-1, 6 * qlsc.nside * qlsc.nside]]

def assert_whole_sphere(fulls:np.ndarray, partials:np.ndarray, qlsc:QLSC, description:str):
	'''
	Assert that a (fulls, partials) pair is q3c's whole-sphere answer.
	'''
	assert fulls.shape == (0,2), f"{description}: expected no fully covered ranges, got {fulls.tolist()}"
	assert partials.tolist() == whole_sphere_partials(qlsc), \
		f"{description}: expected the whole sphere as one partial range, got {partials.tolist()}"

def test_q3c_maxrad_is_exported_and_plausible():
	'''
	Test that the threshold is published by the C extension, so the Python layer can
	reason about it instead of repeating the number.
	'''
	assert isinstance(q3c.Q3C_MAXRAD, float)
	assert 0 < q3c.Q3C_MAXRAD < 180

def test_maxrad_is_below_the_ellipse_limit():
	'''
	Test the relationship the ellipse candidate search depends on.

	The widening loop grows an ellipse only as far as Q3C_MAXRAD, because past that
	point the answer cannot get any wider. That is only a legal ceiling while
	Q3C_MAXRAD is below the largest ellipse the wrapper will accept; if an upstream
	sync raised Q3C_MAXRAD above it, the widening loop would ask for a region the C
	layer refuses and the query would start raising instead of answering.
	'''
	assert q3c.Q3C_MAXRAD < ELLIPSE_MAXIMUM_SEMI_MAJOR_AXIS
	assert ELLIPSE_MAXIMUM_SEMI_MAJOR_AXIS == q3c.ELLIPSE_MAXIMUM_SEMI_MAJOR_AXIS

@pytest.mark.parametrize("depth", [0, 5, 17, 30])
def test_radial_query_falls_back_to_the_whole_sphere(depth):
	'''
	Test that a cone search at or above Q3C_MAXRAD returns the whole sphere.
	'''
	qlsc = QLSC(depth=depth)
	for radius in (q3c.Q3C_MAXRAD, q3c.Q3C_MAXRAD + 0.1, 90., 179., 180.):
		fulls, partials = q3c.radial_query(qlsc._hprm, 123.4, -45.6, radius)
		assert_whole_sphere(fulls, partials, qlsc, f"radial_query(radius={radius}) at depth {depth}")

@pytest.mark.parametrize("depth", [0, 5, 17, 30])
def test_ellipse_query_falls_back_to_the_whole_sphere(depth):
	'''
	Test that an ellipse search at or above Q3C_MAXRAD returns the whole sphere, for
	every semi-major axis the wrapper still accepts.

	This is the one the guard added to the C wrapper depends on: an elliptical region
	between Q3C_MAXRAD and the wrapper's limit is not computed geometrically at all, so
	the fact that its coverage is a superset comes entirely from here.
	'''
	qlsc = QLSC(depth=depth)
	largest = np.nextafter(ELLIPSE_MAXIMUM_SEMI_MAJOR_AXIS, 0.0)
	for semi_major_axis in (q3c.Q3C_MAXRAD, q3c.Q3C_MAXRAD + 0.1, 60., 89.9, largest):
		for axis_ratio in (0.05, 0.5, 1.0):
			fulls, partials = q3c.ellipse_query(qlsc._hprm, 123.4, -45.6,
												semi_major_axis, axis_ratio, 37.)
			assert_whole_sphere(fulls, partials, qlsc,
								f"ellipse_query(semi_major_axis={semi_major_axis}, "
								f"axis_ratio={axis_ratio}) at depth {depth}")

def test_just_below_maxrad_is_computed_geometrically():
	'''
	Test that the fallback is a fallback and not the only path: just below the
	threshold, both queries return real ranges rather than the whole sphere.

	Without this the tests above would still pass if q3c returned the whole sphere for
	every query, which would be a catastrophic (if safe) loss of the whole point of the
	index.
	'''
	just_below = q3c.Q3C_MAXRAD - 0.1

	radial_fulls, radial_partials = q3c.radial_query(q30._hprm, 123.4, -45.6, just_below)
	assert radial_partials.tolist() != whole_sphere_partials(q30)
	assert len(radial_fulls) + len(radial_partials) > 0

	ellipse_fulls, ellipse_partials = q3c.ellipse_query(q30._hprm, 123.4, -45.6, just_below, 0.5, 37.)
	assert ellipse_partials.tolist() != whole_sphere_partials(q30)
	assert len(ellipse_fulls) + len(ellipse_partials) > 0

def test_the_whole_sphere_range_really_covers_every_pixel():
	'''
	Test that the fallback's single range covers the whole scheme, rather than merely
	being a large number that looks right.
	'''
	for depth in (0, 3, 10, 30):
		qlsc = QLSC(depth=depth)
		fulls, partials = q3c.radial_query(qlsc._hprm, 0., 0., 180.)
		low, high = partials[0]
		assert low < 0
		assert high == qlsc.nbins
		# the extremes of the scheme, and a scattering of real positions
		assert low < 0 < qlsc.nbins - 1 < high
		for ra, dec in ((0., 0.), (359.9, 89.9), (180., -89.9), (45., 35.26), (270., -12.)):
			assert low <= qlsc.ang2ipix(ra, dec) < high

# semi_major_axis
refused_sizes = [ELLIPSE_MAXIMUM_SEMI_MAJOR_AXIS, 90.000001, 120., 180.]

@pytest.mark.parametrize("semi_major_axis", refused_sizes)
def test_both_c_level_ellipse_functions_refuse_an_oversized_region(semi_major_axis):
	'''
	Test that the size limit is enforced symmetrically on both exported ellipse
	functions.

	The predicate had the guard and the range query did not, so a caller reaching
	qlsc.q3c.ellipse_query directly got a silent whole-sphere answer for a region that
	is not the one they asked for - safe only because of the fallback above, which is a
	property of the vendored C rather than of this wrapper.
	'''
	with pytest.raises(ValueError):
		q3c.ellipse_query(q30._hprm, 10., 20., semi_major_axis, 0.5, 0.)
	with pytest.raises(ValueError):
		q3c.in_ellipse(ra=1., dec=0., ra_cen=0., dec_cen=0.,
					   semi_major_axis=semi_major_axis, axis_ratio=0.5, position_angle=0.)

def test_the_guard_does_not_block_a_circle():
	'''
	Test that a circle of any radius is still answered: the Python layer routes an axis
	ratio of 1 to the cone search, which has no such limit.
	'''
	from qlsc import QLSCIndex

	index = QLSCIndex(qlsc=QLSC(depth=8))
	index.add_point(ra=110., dec=20., key="far")     # 92.08 degrees from the center
	index.add_point(ra=150., dec=0., key="further")  # 136.04 degrees from the center
	index.add_point(ra=10., dec=20., key="near")     # the center itself

	# a radius past the point where q3c's own conic gives up, and past 90 degrees,
	# where it rejected everything outright
	assert set(index.ellipse_query(ra=10., dec=20., semi_major_axis=120., axis_ratio=1.,
								   return_key=True).key) == {"far", "near"}
	assert set(index.ellipse_query(ra=10., dec=20., semi_major_axis=180., axis_ratio=1.,
								   return_key=True).key) == {"far", "further", "near"}
	with pytest.raises(ValueError):
		index.ellipse_query(ra=10., dec=20., semi_major_axis=120., axis_ratio=0.999)
