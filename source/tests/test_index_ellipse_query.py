
import pytest
import numpy as np

from qlsc import QLSC, QLSCIndex
from qlsc.generate import sunflower_points_on_sphere

# Tests for QLSCIndex.ellipse_query.
#
# The reference answer is a brute-force sweep of every indexed point through
# QLSC.point_in_ellipse, the same exact test the query uses to filter its candidates.
# That shares a predicate with the thing under test, so the circular case is also
# checked against the great-circle distance, an oracle that shares nothing with it -
# see test_circular_ellipse_against_the_distance_oracle and its neighbours below.

sample_points = sunflower_points_on_sphere(n=10_000)

q30 = QLSC(depth=30)

index = QLSCIndex(qlsc=QLSC(depth=10))
index.add_points(points=sample_points)

def as_pairs(result:np.ndarray) -> np.ndarray:
	'''
	Reshape a query result to (n,2).

	These queries already promise (n,2) (see the shape contract tests below); this
	keeps the comparison helpers usable for the brute-force reference arrays too,
	which are built with reshape(-1,2) and can be empty.
	'''
	result = np.asarray(result)
	if result.size == 0:
		return np.zeros((0,2))
	return np.atleast_2d(result)

def sorted_pairs(points:np.ndarray) -> np.ndarray:
	'''
	Sort an (n,2) array of coordinates so two sets can be compared elementwise.
	'''
	points = as_pairs(points)
	return points[np.lexsort((points[:,1], points[:,0]))]

def brute_force_ellipse(ra, dec, semi_major_axis, axis_ratio, position_angle) -> np.ndarray:
	'''
	Return the sample points inside the ellipse, tested one at a time with the exact
	spherical containment test.
	'''
	inside = [p for p in sample_points
			  if q30.point_in_ellipse(ra=p[0], dec=p[1], ra_cen=ra, dec_cen=dec,
									  semi_major_axis=semi_major_axis, axis_ratio=axis_ratio,
									  position_angle=position_angle)]
	return np.array(inside).reshape(-1, 2)

# label, ra, dec, semi_major_axis, axis_ratio, position_angle
ellipse_cases = [
	("mid sky", 10., 20., 5., 0.4, 35.),
	("circle", 200., -45., 6., 1.0, 0.),
	("needle", 200., -45., 8., 0.05, 123.),
	("north pole", 0., 90., 5., 0.5, 20.),
	("south pole", 0., -90., 5., 0.5, 20.),
	("near pole", 123., 89.5, 3., 0.3, 70.),
	("dateline", 0., 0., 4., 0.4, 90.),
	("cube face boundary", 45., 0., 4., 0.6, 45.),
	("tiny", 45., 30., 1e-5, 0.5, 0.),
	("large", 200., -45., 30., 0.6, 10.),
	("larger than the q3c whole-sky threshold", 200., -45., 40., 0.6, 10.),
	("zero size", 200., -45., 0., 1.0, 0.)
]

@pytest.mark.parametrize("label, ra, dec, semi_major_axis, axis_ratio, position_angle",
						 ellipse_cases, ids=[c[0] for c in ellipse_cases])
def test_ellipse_query_matches_brute_force(label, ra, dec, semi_major_axis, axis_ratio, position_angle):
	'''
	Test that an ellipse query returns exactly the indexed points inside the ellipse.
	'''
	matches = sorted_pairs(index.ellipse_query(ra=ra, dec=dec, semi_major_axis=semi_major_axis,
											   axis_ratio=axis_ratio, position_angle=position_angle))
	expected = sorted_pairs(brute_force_ellipse(ra, dec, semi_major_axis, axis_ratio, position_angle))

	assert matches.shape == expected.shape, f"the {label} ellipse query returned a different number of matches than a brute-force calculation"
	np.testing.assert_allclose(matches, expected, err_msg=f"the {label} ellipse query does not agree with a brute-force calculation")

# ra, dec, radius
circle_cases = [
	(10., 20., 3.),
	(0., 0., 1.),
	(0., 90., 7.),
	(0., -90., 7.),
	(45., 0., 5.),
	(359.5, -0.5, 2.),
	(200., -45., 34.9)
]

@pytest.mark.parametrize("ra, dec, radius", circle_cases)
def test_ellipse_with_axis_ratio_one_agrees_with_radial_query(ra, dec, radius):
	'''
	Test that an ellipse degenerating to a circle returns the same points as a cone
	search of the same radius. ellipse_query answers this case BY running the cone
	search, so what this pins is that it keeps doing so - including that the position
	angle, which a circle has no use for, does not change the answer.
	'''
	for position_angle in (0., 42.):
		ellipse_matches = sorted_pairs(index.ellipse_query(ra=ra, dec=dec, semi_major_axis=radius,
														   axis_ratio=1.0, position_angle=position_angle))
		radial_matches = sorted_pairs(index.radial_query(ra=ra, dec=dec, radius=radius))

		assert ellipse_matches.shape == radial_matches.shape, "a circular ellipse query returned a different number of matches than the equivalent radial query"
		np.testing.assert_allclose(ellipse_matches, radial_matches)

# An independent oracle for the circular case. Everything above compares
# ellipse_query with point_in_ellipse, which is the predicate it filters with; these
# compare it with the great-circle distance instead, which shares no code with it.
# Radii at and beyond 90 degrees are the cases that were wrong: Q3C's elliptical
# predicate is a conic in the tangent plane at the center, so it rejected every point
# more than 90 degrees away whatever the radius, and was not monotonic beyond that -
# a 150 degree "circle" held less than a 60 degree one.

def offset_point(ra0:float, dec0:float, separation:float, position_angle:float):
	'''
	Return the point `separation` degrees from (ra0,dec0) along the great circle at
	the given position angle (measured from north through east).
	'''
	d0 = np.radians(dec0)
	s = np.radians(separation)
	pa = np.radians(position_angle)
	dec = np.arcsin(np.sin(d0)*np.cos(s) + np.cos(d0)*np.sin(s)*np.cos(pa))
	dra = np.arctan2(np.sin(pa)*np.sin(s)*np.cos(d0), np.cos(s) - np.sin(d0)*np.sin(dec))
	return (np.degrees(dra) + ra0) % 360, np.degrees(dec)

@pytest.mark.parametrize("radius", [1., 34.9, 35., 60., 89.9, 90., 90.1, 120., 179., 180.])
def test_circular_ellipse_against_the_distance_oracle(radius):
	'''
	Test a circular ellipse query against the angular separation of every indexed
	point, computed directly. This is the check that would have caught the semi-major
	axis being run through Q3C's conic: at radius=120 the conic returned nothing beyond
	90 degrees, so more than half the matches were missing.
	'''
	from qlsc import distance

	matches = sorted_pairs(index.ellipse_query(ra=10., dec=20., semi_major_axis=radius, axis_ratio=1.))
	expected = sorted_pairs(np.array([p for p in sample_points if distance(p[0], p[1], 10., 20.) <= radius]).reshape(-1,2))

	assert matches.shape == expected.shape, f"a circular ellipse of radius {radius} deg did not match the distance oracle"
	np.testing.assert_allclose(matches, expected)

@pytest.mark.parametrize("radius", [0.5, 5., 45., 90., 135.])
def test_circular_ellipse_includes_points_exactly_on_its_edge(radius):
	'''
	Test the boundary convention for a circle: a point at exactly the radius is inside,
	as it is for a cone search. Q3C's conic uses a strict inequality and would have
	excluded it.

	The point is placed due north of a center on the equator, where the separation is a
	difference of declinations and the comparison is exact rather than a round trip
	through a spherical offset.
	'''
	edge_index = QLSCIndex(qlsc=QLSC(depth=20))
	edge_index.add_point(ra=10., dec=radius, key="on-edge")
	edge_index.add_point(ra=10., dec=radius - 1e-6, key="inside")
	edge_index.add_point(ra=10., dec=radius + 1e-6, key="outside")

	found = set(edge_index.ellipse_query(ra=10., dec=0., semi_major_axis=radius,
										 axis_ratio=1., return_key=True).key)
	assert found == {"on-edge", "inside"}
	# and the cone search agrees, point for point
	assert found == set(edge_index.radial_query(ra=10., dec=0., radius=radius, return_key=True).key)

@pytest.mark.parametrize("radius", [0.5, 45., 90., 135.])
def test_circular_ellipse_agrees_with_the_cone_search_off_the_meridian(radius):
	'''
	Test that the circular case matches the cone search for points placed at an
	arbitrary position angle, where a point nominally on the edge may land on either
	side of it by a floating point hair. Whichever side it lands on, the two queries
	must land there together.
	'''
	offset_index = QLSCIndex(qlsc=QLSC(depth=20))
	for i, scale in enumerate((1 - 1e-9, 1.0, 1 + 1e-9)):
		point_ra, point_dec = offset_point(10., 20., radius * scale, 55.)
		offset_index.add_point(ra=point_ra, dec=point_dec, key=f"point-{i}")

	assert (set(offset_index.ellipse_query(ra=10., dec=20., semi_major_axis=radius,
										   axis_ratio=1., return_key=True).key)
			== set(offset_index.radial_query(ra=10., dec=20., radius=radius, return_key=True).key))

def test_point_in_ellipse_circle_against_the_distance_oracle():
	'''
	Test the circular case of the predicate itself against the great-circle distance,
	over random points and radii spanning the whole [0,180] range.
	'''
	from qlsc import distance

	rng = np.random.default_rng(20260830)
	centers = np.column_stack((rng.uniform(0, 360, 40), np.degrees(np.arcsin(rng.uniform(-1, 1, 40)))))
	points = np.column_stack((rng.uniform(0, 360, 40), np.degrees(np.arcsin(rng.uniform(-1, 1, 40)))))

	for (center_ra, center_dec), (point_ra, point_dec) in zip(centers, points):
		separation = distance(point_ra, point_dec, center_ra, center_dec)
		for radius in (0., 1., 30., 89., 90., 91., 120., 179., 180.):
			expected = separation <= radius
			got = q30.point_in_ellipse(ra=point_ra, dec=point_dec, ra_cen=center_ra, dec_cen=center_dec,
									   semi_major_axis=radius, axis_ratio=1., position_angle=17.)
			assert got == expected, (f"a point {separation} degrees from the center was reported "
									 f"{'inside' if got else 'outside'} a circle of radius {radius}")

@pytest.mark.parametrize("depth", [0, 2, 5, 10, 20, 30])
def test_ellipse_query_at_every_depth(depth):
	'''
	Test that the answer does not depend on the depth of the index; the candidate
	ranges are always computed at depth 30 and mapped down.
	'''
	low_depth_index = QLSCIndex(qlsc=QLSC(depth=depth))
	low_depth_index.add_points(points=sample_points)

	matches = sorted_pairs(low_depth_index.ellipse_query(ra=10., dec=20., semi_major_axis=5.,
														 axis_ratio=0.4, position_angle=35.))
	expected = sorted_pairs(brute_force_ellipse(10., 20., 5., 0.4, 35.))

	assert matches.shape == expected.shape, f"the depth={depth} ellipse query returned a different number of matches than a brute-force calculation"
	np.testing.assert_allclose(matches, expected)

def test_ellipse_query_returns_distance():
	'''
	Test that return_distance gives the angular separation of each match from the
	center of the ellipse, in degrees.
	'''
	matches = index.ellipse_query(ra=10., dec=20., semi_major_axis=5., axis_ratio=0.4,
								  position_angle=35., return_distance=True)
	assert set(matches.dtype.names) == {"ra", "dec", "distance"}
	assert len(matches) > 0

	from qlsc import distance
	expected = [distance(ra, dec, 10., 20.) for ra, dec in zip(matches.ra, matches.dec)]
	np.testing.assert_allclose(matches.distance, expected)
	# every match is within the semi-major axis of the center
	assert (matches.distance <= 5.).all()

def test_ellipse_query_returns_keys_and_distance():
	'''
	Test that both flags together give a recarray with all four fields.
	'''
	keyed_index = QLSCIndex(qlsc=QLSC(depth=10))
	keyed_index.add_points(points=sample_points, keys=[f"point-{i}" for i in range(len(sample_points))])

	matches = keyed_index.ellipse_query(ra=10., dec=20., semi_major_axis=5., axis_ratio=0.4,
										position_angle=35., return_key=True, return_distance=True)
	assert set(matches.dtype.names) == {"ra", "dec", "key", "distance"}
	assert len(matches) == len(brute_force_ellipse(10., 20., 5., 0.4, 35.))
	for ra, dec, key in zip(matches.ra, matches.dec, matches.key):
		assert key == f"point-{np.where((sample_points[:,0] == ra) & (sample_points[:,1] == dec))[0][0]}"

def test_ellipse_query_with_no_matches():
	'''
	Test that an ellipse containing none of the indexed points returns an empty result.
	'''
	empty_index = QLSCIndex(qlsc=QLSC(depth=10))
	empty_index.add_point(ra=180., dec=0.)

	matches = empty_index.ellipse_query(ra=10., dec=20., semi_major_axis=1., axis_ratio=0.5)
	assert as_pairs(matches).shape == (0,2)

def test_zero_size_ellipse_finds_equivalent_coordinate_representations():
	'''
	Test that a zero-size ellipse behaves as a zero-radius cone search does: it matches
	a point stored as an equivalent representation of its center (ra=360 rather than
	ra=0, or any ra at a pole), where the exact test is defeated by floating point
	residuals and q3c's candidate ranges do not reach the neighbouring pixel.
	'''
	seam_index = QLSCIndex(qlsc=QLSC(depth=30))
	seam_index.add_point(ra=360., dec=20., key="seam")
	seam_index.add_point(ra=123., dec=90., key="pole")

	assert set(seam_index.ellipse_query(ra=0., dec=20., semi_major_axis=0., axis_ratio=1.,
										return_key=True).key) == {"seam"}
	assert set(seam_index.ellipse_query(ra=0., dec=90., semi_major_axis=0., axis_ratio=1.,
										return_key=True).key) == {"pole"}
	# and the reported separation of such a match is exactly zero
	assert seam_index.ellipse_query(ra=0., dec=20., semi_major_axis=0., axis_ratio=1.,
									return_distance=True).distance[0] == 0.

def test_ellipse_query_requires_its_parameters():
	'''
	Test that omitting any required parameter raises rather than being silently defaulted.
	'''
	full = dict(ra=10., dec=20., semi_major_axis=1., axis_ratio=0.5)
	for name in full:
		partial = dict(full)
		partial[name] = None
		with pytest.raises(ValueError):
			index.ellipse_query(**partial)

# parameter, value
invalid_parameter_cases = [
	("axis_ratio", 0.),
	("axis_ratio", -0.5),
	("axis_ratio", 1.5),
	("axis_ratio", np.nan),
	("semi_major_axis", -1.),
	("semi_major_axis", 181.),
	("semi_major_axis", np.inf),
	("position_angle", np.nan),
	("ra", np.nan),
	("dec", 91.),
	("dec", np.nan)
]

@pytest.mark.parametrize("parameter, value", invalid_parameter_cases)
def test_invalid_parameters_raise(parameter, value):
	'''
	Test that out of range and non-finite parameters raise ValueError.
	'''
	arguments = dict(ra=10., dec=20., semi_major_axis=1., axis_ratio=0.5, position_angle=0.)
	arguments[parameter] = value
	with pytest.raises(ValueError):
		index.ellipse_query(**arguments)

def test_astropy_quantity_parameters():
	'''
	Test that the semi-major axis and position angle accept astropy Quantities, as
	radial_query's radius does.
	'''
	astropy_units = pytest.importorskip("astropy.units")
	quantities = index.ellipse_query(ra=10., dec=20., semi_major_axis=300 * astropy_units.arcmin,
									 axis_ratio=0.4, position_angle=np.pi/4 * astropy_units.rad)
	degrees = index.ellipse_query(ra=10., dec=20., semi_major_axis=5., axis_ratio=0.4, position_angle=45.)
	np.testing.assert_allclose(sorted_pairs(quantities), sorted_pairs(degrees))

def test_ellipse_query_result_is_always_two_dimensional():
	'''
	Test the return shape contract: (n,2) for any number of matches, including one.

	radial_query squeezes its result, so a single match comes back with shape (2,) and
	``matches[:,0]`` fails for exactly that case. That wart is kept for compatibility
	and queued for a 2.0 release; this query must not inherit it, and neither must the
	circular case, which is answered by handing the query to radial_query.
	'''
	one_point = QLSCIndex(qlsc=QLSC(depth=10))
	one_point.add_point(ra=10., dec=20.)

	for axis_ratio in (0.5, 1.0):
		matches = one_point.ellipse_query(ra=10., dec=20., semi_major_axis=1., axis_ratio=axis_ratio)
		assert matches.shape == (1,2)
		assert matches[:,0][0] == pytest.approx(10.)

	empty = QLSCIndex(qlsc=QLSC(depth=10))
	empty.add_point(ra=180., dec=0.)
	assert empty.ellipse_query(ra=10., dec=20., semi_major_axis=1., axis_ratio=0.5).shape == (0,2)
	assert empty.ellipse_query(ra=10., dec=20., semi_major_axis=1., axis_ratio=1.0).shape == (0,2)

	many = index.ellipse_query(ra=10., dec=20., semi_major_axis=5., axis_ratio=0.4)
	assert many.ndim == 2 and many.shape[1] == 2

def test_radial_query_keeps_its_historical_shape():
	'''
	Test that the shape fix above was NOT applied to radial_query; changing it would
	break existing callers, and it is queued for 2.0 instead.
	'''
	one_point = QLSCIndex(qlsc=QLSC(depth=10))
	one_point.add_point(ra=10., dec=20.)
	assert np.asarray(one_point.radial_query(ra=10., dec=20., radius=1.)).shape == (2,)

def test_ellipse_query_on_a_closed_index_raises():
	'''
	Test that querying a closed index raises rather than reopening it silently.
	'''
	closed_index = QLSCIndex(qlsc=QLSC(depth=10))
	closed_index.add_point(ra=10., dec=20.)
	closed_index.close()
	with pytest.raises(Exception):
		closed_index.ellipse_query(ra=10., dec=20., semi_major_axis=1., axis_ratio=0.5)
