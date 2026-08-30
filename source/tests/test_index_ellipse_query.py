
import pytest
import numpy as np

from qlsc import QLSC, QLSCIndex
from qlsc.generate import sunflower_points_on_sphere

# Tests for QLSCIndex.ellipse_query.
#
# The reference answer is a brute-force sweep of every indexed point through
# QLSC.point_in_ellipse, the same exact test the query uses to filter its
# candidates; an ellipse with an axis ratio of 1 is additionally checked against
# radial_query, which reaches the same points by an entirely different route.

sample_points = sunflower_points_on_sphere(n=10_000)

q30 = QLSC(depth=30)

index = QLSCIndex(qlsc=QLSC(depth=10))
index.add_points(points=sample_points)

def as_pairs(result:np.ndarray) -> np.ndarray:
	'''
	Reshape a query result to (n,2). A single match is returned with shape (2,) and
	no matches with shape (0,2); this is the documented (and awkward) return shape.
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
	search of the same radius. The two share no filtering code: radial_query compares
	sin^2(d/2) against the radius, ellipse_query evaluates q3c's conic.
	'''
	for position_angle in (0., 42.):
		ellipse_matches = sorted_pairs(index.ellipse_query(ra=ra, dec=dec, semi_major_axis=radius,
														   axis_ratio=1.0, position_angle=position_angle))
		radial_matches = sorted_pairs(index.radial_query(ra=ra, dec=dec, radius=radius))

		assert ellipse_matches.shape == radial_matches.shape, "a circular ellipse query returned a different number of matches than the equivalent radial query"
		np.testing.assert_allclose(ellipse_matches, radial_matches)

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

def test_ellipse_query_on_a_closed_index_raises():
	'''
	Test that querying a closed index raises rather than reopening it silently.
	'''
	closed_index = QLSCIndex(qlsc=QLSC(depth=10))
	closed_index.add_point(ra=10., dec=20.)
	closed_index.close()
	with pytest.raises(Exception):
		closed_index.ellipse_query(ra=10., dec=20., semi_major_axis=1., axis_ratio=0.5)
