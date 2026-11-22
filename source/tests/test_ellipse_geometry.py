
import pytest
import numpy as np

from qlsc import QLSC, QLSCIndex, q3c
from qlsc.qlsc_base import ELLIPSE_MAXIMUM_SEMI_MAJOR_AXIS

# What Q3C's ellipse actually IS, pinned as a test rather than left as folklore.
#
# q3c_in_ellipse evaluates a conic in the gnomonic (tangent plane) projection about
# the center of the ellipse. In that plane the region has semi-axes tan(a) and
# tan(a)*axis_ratio, so:
#
#   * along the major axis the edge is at exactly `a` degrees, at any size;
#   * across the minor axis the edge is at arctan(tan(a)*axis_ratio), which is close
#     to a*axis_ratio for a small ellipse and much wider for a large one (at a=80 and
#     axis_ratio=0.2 the semi-minor axis is 48.6 degrees, not 16);
#   * only the near hemisphere has a tangent plane, so the whole construction stops
#     meaning anything at 90 degrees. Past it q3c silently answers about a different
#     region, which is why ELLIPSE_MAXIMUM_SEMI_MAJOR_AXIS exists.
#
# This was read out of q3cube.c and then verified numerically; these tests are that
# verification, kept so that an upstream sync cannot change the shape unnoticed.

q30 = QLSC(depth=30)

def offset_point(ra0:float, dec0:float, separation:float, position_angle:float):
	'''
	Return the point `separation` degrees from (ra0,dec0) along the great circle at
	the given position angle (measured from north through east).

	Computed with unit vectors and arctan2 rather than arcsin: arcsin cannot resolve a
	declination within ~1e-6 degrees of a pole (its argument rounds to 1), which put
	every probe point that close to a pole on the same few positions and made an
	ellipse whose edge passes through the pole look a microdegree too small.
	'''
	a0, d0 = np.radians(ra0), np.radians(dec0)
	s = np.radians(separation)
	pa = np.radians(position_angle)
	center = np.array([np.cos(d0)*np.cos(a0), np.cos(d0)*np.sin(a0), np.sin(d0)])
	north = np.array([-np.sin(d0)*np.cos(a0), -np.sin(d0)*np.sin(a0), np.cos(d0)])
	east = np.array([-np.sin(a0), np.cos(a0), 0.])
	x, y, z = np.cos(s)*center + np.sin(s)*(np.cos(pa)*north + np.sin(pa)*east)
	return np.degrees(np.arctan2(y, x)) % 360, np.degrees(np.arctan2(z, np.hypot(x, y)))

def edge_separation(ra0:float, dec0:float, semi_major_axis:float, axis_ratio:float,
					position_angle:float, direction:float) -> float:
	'''
	Find, by bisection, how far the edge of the ellipse is from its center along the
	great circle at position angle `direction`.
	'''
	low, high = 1e-9, 89.9999999
	for _ in range(80):
		middle = (low + high) / 2
		ra, dec = offset_point(ra0, dec0, middle, direction)
		if q30.point_in_ellipse(ra=ra, dec=dec, ra_cen=ra0, dec_cen=dec0,
								semi_major_axis=semi_major_axis, axis_ratio=axis_ratio,
								position_angle=position_angle):
			low = middle
		else:
			high = middle
	return low

def tangent_plane_edge(semi_major_axis:float, axis_ratio:float, azimuth:float) -> float:
	'''
	The separation of the ellipse's edge from its center at `azimuth` degrees from the
	major axis, if the ellipse is the tangent plane conic with semi-axes tan(a) and
	tan(a)*axis_ratio.
	'''
	ta = np.tan(np.radians(semi_major_axis))
	phi = np.radians(azimuth)
	return np.degrees(np.arctan(1.0 / np.hypot(np.cos(phi)/ta, np.sin(phi)/(ta*axis_ratio))))

# center_ra, center_dec
centers = [(0., 0.), (123.4, -67.8), (310., 85.), (45., 20.)]

@pytest.mark.parametrize("semi_major_axis", [0.5, 5., 20., 34.9, 50., 75., 89.])
@pytest.mark.parametrize("axis_ratio", [0.05, 0.3, 0.7, 0.999])
def test_ellipse_is_the_tangent_plane_conic(semi_major_axis, axis_ratio):
	'''
	Test that the edge of the ellipse lies where the gnomonic tangent plane conic puts
	it, in every direction, at every center and position angle.
	'''
	for center_ra, center_dec in centers:
		for position_angle in (0., 37., 90., 210.):
			for azimuth in (0., 23., 45., 67., 90.):
				measured = edge_separation(center_ra, center_dec, semi_major_axis, axis_ratio,
										   position_angle, position_angle + azimuth)
				expected = tangent_plane_edge(semi_major_axis, axis_ratio, azimuth)
				assert measured == pytest.approx(expected, abs=1e-8), (
					f"the ellipse edge at azimuth {azimuth} from the major axis is {measured} degrees "
					f"from the center, not the {expected} the tangent plane conic predicts")

@pytest.mark.parametrize("semi_major_axis", [1., 10., 45., 89.])
def test_major_axis_is_exactly_the_semi_major_axis(semi_major_axis):
	'''
	Test that the ellipse reaches exactly its semi-major axis along the major axis,
	whatever the axis ratio; this is the one dimension that is not stretched.
	'''
	for axis_ratio in (0.1, 0.5, 0.9):
		measured = edge_separation(45., 20., semi_major_axis, axis_ratio, 30., 30.)
		assert measured == pytest.approx(semi_major_axis, abs=1e-8)

def test_the_minor_axis_is_stretched_for_a_large_ellipse():
	'''
	Test the consequence that makes the 90 degree limit necessary: the minor axis is
	arctan(tan(a)*ratio), not a*ratio, so as the semi-major axis approaches 90 degrees
	the "ellipse" swells towards a hemisphere whatever its axis ratio.
	'''
	small = edge_separation(0., 0., 1., 0.2, 0., 90.)
	large = edge_separation(0., 0., 80., 0.2, 0., 90.)
	nearly_hemispherical = edge_separation(0., 0., 89.9, 0.2, 0., 90.)

	assert small == pytest.approx(0.2, abs=1e-3)			# a small ellipse behaves linearly
	assert large == pytest.approx(48.599, abs=1e-3)			# arctan(tan(80)*0.2)
	assert nearly_hemispherical > 89.				# a 1:5 "ellipse" that is nearly a circle

# semi_major_axis
too_large_semi_major_axes = [ELLIPSE_MAXIMUM_SEMI_MAJOR_AXIS, 90.000001, 120., 150., 180.]

@pytest.mark.parametrize("semi_major_axis", too_large_semi_major_axes)
def test_elliptical_regions_past_the_hemisphere_are_refused(semi_major_axis):
	'''
	Test that a genuinely elliptical region larger than the tangent plane can describe
	raises, at every layer, rather than answering about some other region.

	Without the limit, q3c_in_ellipse rejects every point more than 90 degrees from the
	center whatever the semi-major axis, and the sin^2 term it tests is not monotonic
	past 90 degrees, so a 150 degree "ellipse" is SMALLER than a 60 degree one.
	'''
	with pytest.raises(ValueError):
		q30.point_in_ellipse(ra=1., dec=0., ra_cen=0., dec_cen=0.,
							 semi_major_axis=semi_major_axis, axis_ratio=0.5)

	with pytest.raises(ValueError):
		q3c.in_ellipse(ra=1., dec=0., ra_cen=0., dec_cen=0.,
					   semi_major_axis=semi_major_axis, axis_ratio=0.5, position_angle=0.)

	index = QLSCIndex(qlsc=QLSC(depth=10))
	index.add_point(ra=1., dec=0.)
	with pytest.raises(ValueError):
		index.ellipse_query(ra=0., dec=0., semi_major_axis=semi_major_axis, axis_ratio=0.5)

def test_just_below_the_limit_is_accepted():
	'''
	Test that the limit is not over-tight: the geometry is exact right up to it.
	'''
	assert q30.point_in_ellipse(ra=1., dec=0., ra_cen=0., dec_cen=0.,
								semi_major_axis=89.999, axis_ratio=0.5) is True

def test_circles_are_not_limited_to_the_hemisphere():
	'''
	Test that an axis ratio of exactly 1 is exempt: a circle is a distance test, needs
	no tangent plane, and is accepted over the full [0,180] range a cone search takes.
	'''
	for semi_major_axis in (90., 120., 179., 180.):
		assert q30.point_in_ellipse(ra=90., dec=0., ra_cen=0., dec_cen=0.,
									semi_major_axis=semi_major_axis, axis_ratio=1.) is True
	with pytest.raises(ValueError):
		q30.point_in_ellipse(ra=90., dec=0., ra_cen=0., dec_cen=0.,
							 semi_major_axis=180.1, axis_ratio=1.)
