
import pytest
import numpy as np

from qlsc import QLSC, QLSCIndex, q3c
from qlsc.generate import sunflower_points_on_sphere

# Tests for QLSCIndex.polygon_query.
#
# The reference answer everywhere is a brute-force sweep of every indexed point
# through QLSC.point_in_polygon: the query is only allowed to be a faster way of
# computing the same set, never a different one.

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

def brute_force_polygon(polygon:np.ndarray) -> np.ndarray:
	'''
	Return the sample points inside the polygon, tested one at a time with the exact
	spherical containment test.
	'''
	inside = [p for p in sample_points if q30.point_in_polygon(ra=p[0], dec=p[1], polygon=polygon)]
	return np.array(inside).reshape(-1, 2)

# The polygons are given as [ra,dec] vertex pairs with no repeated endpoint.
small_square = np.array([[9.,19.],[11.,19.],[11.,21.],[9.,21.]])
face_boundary = np.array([[40.,-5.],[50.,-5.],[50.,5.],[40.,5.]])
dateline = np.array([[355.,-5.],[5.,-5.],[5.,5.],[355.,5.]])
north_pole = np.array([[0.,85.],[90.,85.],[180.,85.],[270.,85.]])
south_pole = np.array([[0.,-85.],[90.,-85.],[180.,-85.],[270.,-85.]])
triangle = np.array([[200.,-40.],[210.,-40.],[205.,-32.]])
twenty_degrees = np.array([[100.,-10.],[120.,-10.],[120.,10.],[100.,10.]])
tiny = np.array([[10.,20.],[10.0001,20.],[10.0001,20.0001],[10.,20.0001]])
concave_chevron = np.array([[0.,0.],[10.,0.],[10.,10.],[5.,2.],[0.,10.]])
concave_comb = np.array([[0.,0.],[10.,0.],[10.,10.],[8.,1.],[6.,10.],[4.,1.],[2.,10.],[0.,10.]])

polygon_cases = [
	("small square", small_square),
	("cube face boundary", face_boundary),
	("dateline crossing", dateline),
	("north pole", north_pole),
	("south pole", south_pole),
	("triangle", triangle),
	("twenty degrees across", twenty_degrees),
	("tiny", tiny),
	("concave chevron", concave_chevron),
	("concave comb", concave_comb)
]

@pytest.mark.parametrize("label, polygon", polygon_cases, ids=[c[0] for c in polygon_cases])
def test_polygon_query_matches_brute_force(label, polygon):
	'''
	Test that a polygon query returns exactly the indexed points inside the polygon.
	'''
	matches = sorted_pairs(index.polygon_query(polygon=polygon))
	expected = sorted_pairs(brute_force_polygon(polygon))

	assert matches.shape == expected.shape, f"the {label} polygon query returned a different number of matches than a brute-force calculation"
	np.testing.assert_allclose(matches, expected, err_msg=f"the {label} polygon query does not agree with a brute-force calculation")

@pytest.mark.parametrize("depth", [0, 2, 5, 10, 20, 30])
def test_polygon_query_at_every_depth(depth):
	'''
	Test that the answer does not depend on the depth of the index. The candidate
	ranges are always computed at depth 30 and mapped down, so a low-depth index
	simply has coarser candidate pixels and more of them to filter.
	'''
	low_depth_index = QLSCIndex(qlsc=QLSC(depth=depth))
	low_depth_index.add_points(points=sample_points)

	matches = sorted_pairs(low_depth_index.polygon_query(polygon=small_square))
	expected = sorted_pairs(brute_force_polygon(small_square))

	assert matches.shape == expected.shape, f"the depth={depth} polygon query returned a different number of matches than a brute-force calculation"
	np.testing.assert_allclose(matches, expected)

def test_polygon_query_returns_keys():
	'''
	Test that return_key produces a recarray carrying the keys the points were added with.
	'''
	keyed_index = QLSCIndex(qlsc=QLSC(depth=10))
	keyed_index.add_points(points=sample_points, keys=[f"point-{i}" for i in range(len(sample_points))])

	matches = keyed_index.polygon_query(polygon=small_square, return_key=True)
	expected = brute_force_polygon(small_square)

	assert len(matches) == len(expected)
	assert set(matches.dtype.names) == {"ra", "dec", "key"}
	for ra, dec, key in zip(matches.ra, matches.dec, matches.key):
		assert key == f"point-{np.where((sample_points[:,0] == ra) & (sample_points[:,1] == dec))[0][0]}"

@pytest.mark.parametrize("label, polygon", polygon_cases, ids=[c[0] for c in polygon_cases])
def test_batched_filter_agrees_with_the_scalar_one(label, polygon):
	'''
	Test that the batched containment test the query filters with gives exactly the
	same answers as the scalar one, point for point.

	polygon_query uses check_sphere_points_in_poly, which projects the polygon onto the
	cube faces once and reuses that projection (q3c's invocation=1 path) instead of
	reprojecting every vertex for every candidate. That is a large speed-up and a place
	where a stale or mis-shared projection buffer would show up as wrong answers.
	'''
	poly_ra = np.ascontiguousarray(np.asarray(polygon, dtype=np.double)[:,0])
	poly_dec = np.ascontiguousarray(np.asarray(polygon, dtype=np.double)[:,1])

	ra = np.ascontiguousarray(sample_points[:,0])
	dec = np.ascontiguousarray(sample_points[:,1])

	batched = q3c.check_sphere_points_in_poly(q30._hprm, ra, dec, poly_ra, poly_dec)
	scalar = np.array([q3c.check_sphere_point_in_poly(q30._hprm, r, d, poly_ra, poly_dec)
					   for r, d in zip(ra, dec)])

	assert batched.dtype == np.bool_
	assert batched.shape == (len(sample_points),)
	np.testing.assert_array_equal(batched, scalar,
								  err_msg=f"the batched filter disagrees with the scalar one for the {label} polygon")

def test_batched_filter_edge_cases():
	'''
	Test the batched filter's own boundaries: an empty point array still reports a
	polygon that cannot be projected, mismatched lengths and non-finite points raise.
	'''
	poly_ra = np.ascontiguousarray(small_square[:,0])
	poly_dec = np.ascontiguousarray(small_square[:,1])
	empty = np.array([], dtype=np.double)

	assert q3c.check_sphere_points_in_poly(q30._hprm, empty, empty, poly_ra, poly_dec).shape == (0,)

	huge = np.array([[0.,-30.],[60.,-30.],[60.,30.],[0.,30.]])
	with pytest.raises(ValueError):
		# the too-large polygon must be reported even with no points to test
		q3c.check_sphere_points_in_poly(q30._hprm, empty, empty,
										np.ascontiguousarray(huge[:,0]), np.ascontiguousarray(huge[:,1]))

	with pytest.raises(ValueError):
		q3c.check_sphere_points_in_poly(q30._hprm, np.array([1., 2.]), np.array([1.]), poly_ra, poly_dec)

	with pytest.raises(ValueError):
		q3c.check_sphere_points_in_poly(q30._hprm, np.array([np.nan]), np.array([0.]), poly_ra, poly_dec)

def test_polygon_query_result_is_always_two_dimensional():
	'''
	Test the return shape contract: (n,2) for any number of matches, including one.

	radial_query squeezes its result, so a single match comes back with shape (2,) and
	``matches[:,0]`` fails for exactly that case. That wart is kept for compatibility
	and queued for a 2.0 release; this query must not inherit it.
	'''
	one_point = QLSCIndex(qlsc=QLSC(depth=10))
	one_point.add_point(ra=10., dec=20.)

	matches = one_point.polygon_query(polygon=small_square)
	assert matches.shape == (1,2)
	assert matches[:,0][0] == pytest.approx(10.)

	many = index.polygon_query(polygon=twenty_degrees)
	assert many.ndim == 2 and many.shape[1] == 2

def test_polygon_query_with_no_matches():
	'''
	Test that a polygon containing none of the indexed points returns an empty result
	rather than raising or returning a scalar.
	'''
	empty_index = QLSCIndex(qlsc=QLSC(depth=10))
	empty_index.add_point(ra=180., dec=0.)

	matches = empty_index.polygon_query(polygon=small_square)
	assert as_pairs(matches).shape == (0,2)

def test_polygon_query_finds_equivalent_coordinate_representations():
	'''
	Test that points stored in an equivalent representation of their position are
	found: ra=360 rather than ra=0 at the seam, and any ra at all at a pole.
	'''
	seam_index = QLSCIndex(qlsc=QLSC(depth=30))
	seam_index.add_point(ra=360., dec=0.5, key="seam")
	seam_index.add_point(ra=0.5, dec=0.5, key="plain")
	seam_index.add_point(ra=180., dec=0.5, key="elsewhere")

	matches = seam_index.polygon_query(polygon=dateline, return_key=True)
	assert set(matches.key) == {"seam", "plain"}

	pole_index = QLSCIndex(qlsc=QLSC(depth=30))
	pole_index.add_point(ra=123., dec=90., key="pole")
	assert set(pole_index.polygon_query(polygon=north_pole, return_key=True).key) == {"pole"}

def test_concave_polygon_coarse_coverage_is_generous_but_filtered():
	'''
	Test the concave polygon case explicitly.

	Q3C calls a pixel of the coarse search fully covered when all four of its corners
	are inside the polygon (q3c_poly_cover_check), which for a concave polygon is not
	enough: area in a notch can sit inside such a pixel. Upstream issue #55 records
	this as over-coverage. This test first shows the over-coverage really happens for
	the comb polygon - there are positions outside the polygon whose depth-30 pixel
	lies in a range Q3C reported as FULLY covered - and then shows it makes no
	difference to the answer, because every candidate is put through the exact
	containment test before it is returned.

	Over-coverage is safe; under-coverage would not be. The second half of this test
	is what guards that: any position inside the polygon that the coarse ranges missed
	would be a point the query loses.
	'''
	poly_ra = np.ascontiguousarray(concave_comb[:,0])
	poly_dec = np.ascontiguousarray(concave_comb[:,1])
	fulls, partials = q3c.poly_query(q30._hprm, poly_ra, poly_dec)

	# a grid over the polygon's bounding box, in and out of the comb's notches
	grid_ra, grid_dec = np.meshgrid(np.linspace(-0.5, 10.5, 120), np.linspace(-0.5, 10.5, 120))
	grid_ra = grid_ra.ravel()
	grid_dec = grid_dec.ravel()
	ipix = q30.ang2ipix(grid_ra, grid_dec)
	inside = np.array([q30.point_in_polygon(ra=r, dec=d, polygon=concave_comb) for r,d in zip(grid_ra, grid_dec)])

	def covered_by(ranges):
		hit = np.zeros(len(ipix), dtype=bool)
		for low, high in ranges:
			hit |= (ipix >= low) & (ipix < high)
		return hit

	in_fulls = covered_by(fulls)
	in_any = in_fulls | covered_by(partials)

	assert (~inside & in_fulls).sum() > 0, "expected the concave polygon to be over-covered by q3c's 'fully covered' ranges"
	assert (inside & ~in_any).sum() == 0, "the coarse ranges missed part of the polygon; this would silently lose points"

	# and the query itself still agrees exactly with the brute-force answer
	matches = sorted_pairs(index.polygon_query(polygon=concave_comb))
	expected = sorted_pairs(brute_force_polygon(concave_comb))
	assert matches.shape == expected.shape
	np.testing.assert_allclose(matches, expected)

def test_polygon_too_large_raises():
	'''
	Test that a polygon too large to project onto a cube face raises ValueError rather
	than returning a wrong (partial) answer.
	'''
	# 60 degrees across; point_in_polygon accepts this one, the query cannot
	huge = np.array([[0.,-30.],[60.,-30.],[60.,30.],[0.,30.]])
	with pytest.raises(ValueError):
		index.polygon_query(polygon=huge)

def test_polygon_query_requires_a_polygon():
	'''
	Test that omitting the polygon raises.
	'''
	with pytest.raises(ValueError):
		index.polygon_query()

# polygon, reason
invalid_polygon_cases = [
	(np.array([[0.,0.],[1.,0.]]), "only two vertices"),
	(np.array([[0.,0.,0.],[1.,0.,0.],[1.,1.,0.]]), "vertices are not [ra,dec] pairs"),
	(np.array([0.,1.,2.]), "not a 2D array"),
	(np.array([[0.,0.],[1.,0.],[np.nan,1.]]), "a non-finite vertex"),
	(np.array([[0.,0.],[1.,0.],[1.,91.]]), "a dec outside [-90,90]"),
	(np.zeros((101,2)), "more vertices than Q3C allows")
]

@pytest.mark.parametrize("polygon, reason", invalid_polygon_cases, ids=[c[1] for c in invalid_polygon_cases])
def test_invalid_polygons_raise(polygon, reason):
	'''
	Test that a polygon that cannot be interpreted raises ValueError.
	'''
	with pytest.raises(ValueError):
		index.polygon_query(polygon=polygon)

def test_polygon_query_accepts_a_list_of_pairs():
	'''
	Test that the polygon can be given in any form NumPy will make an (n,2) array of,
	as point_in_polygon accepts.
	'''
	from_list = index.polygon_query(polygon=[[9.,19.],[11.,19.],[11.,21.],[9.,21.]])
	from_array = index.polygon_query(polygon=small_square)
	np.testing.assert_allclose(sorted_pairs(from_list), sorted_pairs(from_array))

def test_polygon_query_on_a_closed_index_raises():
	'''
	Test that querying a closed index raises rather than reopening it silently.
	'''
	closed_index = QLSCIndex(qlsc=QLSC(depth=10))
	closed_index.add_point(ra=10., dec=20.)
	closed_index.close()
	with pytest.raises(Exception):
		closed_index.polygon_query(polygon=small_square)
