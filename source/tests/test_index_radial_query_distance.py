
import pytest
import numpy as np

from qlsc import QLSC, QLSCIndex, distance
from qlsc.generate import sunflower_points_on_sphere

# QLSCIndex.radial_query can return the angular separation of each match from the
# query center. The cone filter already computes sin^2(d/2) for every candidate,
# so the value is free; the flag exists only to keep the default return shape
# backward compatible.


def _index(points, keys=None, depth=30):
	q = QLSC(depth=depth)
	index = QLSCIndex(qlsc=q)
	points = np.asarray(points, dtype=np.double)
	index.add_points(ra=points[:,0], dec=points[:,1], keys=keys)
	return index


def test_default_return_is_unchanged():
	'''
	Test that leaving both flags off returns the historical (n,2) array.
	'''
	index = _index([[10., 30.], [10.5, 30.], [11., 30.], [200., -10.]])
	matches = index.radial_query(10., 30., 1.0)

	assert isinstance(matches, np.ndarray)
	assert matches.dtype.names is None, "the default return must not become a recarray"
	assert matches.shape == (3, 2)


def test_distance_field_added():
	'''
	Test that return_distance produces a recarray with a "distance" field.
	'''
	index = _index([[10., 30.], [10.5, 30.], [11., 30.]])
	matches = index.radial_query(10., 30., 1.0, return_distance=True)

	assert matches.dtype.names == ("ra", "dec", "distance")
	assert len(matches) == 3


def test_distance_values_match_the_distance_function():
	'''
	Test that the reported separations agree with the public qlsc.distance().
	'''
	points = sunflower_points_on_sphere(n=5000)
	index = _index(points)

	for center_ra, center_dec in [(0., 0.), (10., 30.), (359.9, -0.1), (123.4, 89.9), (250., -89.99), (45., 35.26)]:
		matches = index.radial_query(center_ra, center_dec, 15.0, return_distance=True)
		assert len(matches) > 0, f"no matches at ({center_ra},{center_dec}); the test would be vacuous"
		for row in matches:
			expected = distance(center_ra, center_dec, row["ra"], row["dec"])
			assert row["distance"] == pytest.approx(expected, abs=1e-9)


def test_distances_are_within_the_search_radius():
	'''
	Test that every reported separation is inside the requested radius.
	'''
	points = sunflower_points_on_sphere(n=20000)
	index = _index(points)

	for radius in [0.5, 5.0, 45.0, 100.0]:
		matches = index.radial_query(77., -12., radius, return_distance=True)
		assert len(matches) > 0
		assert (matches["distance"] <= radius + 1e-9).all()
		assert (matches["distance"] >= 0).all()


def test_distance_with_keys():
	'''
	Test that return_key and return_distance combine, and that the key values
	are unaffected by the distance calculation.
	'''
	index = _index([[10., 30.], [10.5, 30.], [11., 30.]], keys=["a", "b", "c"])

	matches = index.radial_query(10., 30., 1.0, return_key=True, return_distance=True)
	assert matches.dtype.names == ("ra", "dec", "key", "distance")
	assert sorted(matches["key"].tolist()) == ["a", "b", "c"]

	# the key-only form is unchanged
	keys_only = index.radial_query(10., 30., 1.0, return_key=True)
	assert keys_only.dtype.names == ("ra", "dec", "key")

	by_key = {row["key"]: row["distance"] for row in matches}
	assert by_key["a"] == 0.0
	assert by_key["b"] < by_key["c"]


def test_distance_of_an_exact_match_is_zero():
	'''
	Test that a point coincident with the query center reports exactly zero.
	'''
	index = _index([[10., 30.], [200., -10.]])
	matches = index.radial_query(10., 30., 1.0, return_distance=True)
	assert matches["distance"][0] == 0.0


@pytest.mark.parametrize("stored_ra, stored_dec, query_ra, query_dec", [
	(360., 30., 0., 30.),    # equivalent ra representations
	(0., 30., 360., 30.),
	(17., 90., 250., 90.),   # any ra at the north pole is the same point
	(17., -90., 250., -90.)
])
def test_equivalent_representations_report_zero(stored_ra, stored_dec, query_ra, query_dec):
	'''
	Test that a point stored as an equivalent representation of the query center
	reports a separation of exactly zero rather than a floating point residual.

	These matches bypass the sindist filter entirely (they are matched by
	canonical position), so the reported distance is set, not computed.
	'''
	index = _index([[stored_ra, stored_dec]])
	matches = index.radial_query(query_ra, query_dec, 0.0, return_distance=True)
	assert len(matches) == 1
	assert matches["distance"][0] == 0.0


def test_zero_radius_query_with_distance():
	'''
	Test that a zero-radius query still works with return_distance (the sindist
	value is never computed on that path).
	'''
	index = _index([[10., 30.], [10.0001, 30.]])
	matches = index.radial_query(10., 30., 0.0, return_distance=True)
	assert len(matches) == 1
	assert matches["distance"][0] == 0.0


def test_antipodal_distance():
	'''
	Test the extreme of the domain: a 180 degree search reporting an antipodal
	point, where sin^2(d/2) rounds to just over 1.
	'''
	index = _index([[10., 30.], [190., -30.]])
	matches = index.radial_query(10., 30., 180.0, return_distance=True)

	assert len(matches) == 2
	assert max(matches["distance"]) == pytest.approx(180.0, abs=1e-9)
	assert np.isfinite(matches["distance"]).all()


def test_empty_result_shape():
	'''
	Test that a query with no matches returns an empty recarray with the fields
	still present, rather than raising or returning something unindexable.
	'''
	index = _index([[10., 30.]])

	matches = index.radial_query(200., -40., 1.0, return_distance=True)
	assert len(matches) == 0
	assert matches.dtype.names == ("ra", "dec", "distance")

	matches = index.radial_query(200., -40., 1.0, return_key=True, return_distance=True)
	assert len(matches) == 0
	assert matches.dtype.names == ("ra", "dec", "key", "distance")


@pytest.mark.parametrize("depth", [1, 5, 15, 30])
def test_distance_at_every_index_depth(depth):
	'''
	Test that the separation is correct regardless of the index depth (the
	candidate ranges are always computed at depth 30 and mapped down).
	'''
	points = sunflower_points_on_sphere(n=2000)
	index = _index(points, depth=depth)

	matches = index.radial_query(45., 35., 10.0, return_distance=True)
	assert len(matches) > 0
	for row in matches:
		assert row["distance"] == pytest.approx(distance(45., 35., row["ra"], row["dec"]), abs=1e-9)


def test_distance_flag_does_not_change_which_points_match():
	'''
	Test that asking for the distance does not change the match set - the
	restructured filter must behave exactly as it did before.
	'''
	points = sunflower_points_on_sphere(n=10000)
	index = _index(points)

	for center_ra, center_dec, radius in [(0., 0., 3.), (180., 45., 0.5), (270., -89., 12.), (10., 30., 0.)]:
		plain = index.radial_query(center_ra, center_dec, radius)
		with_distance = index.radial_query(center_ra, center_dec, radius, return_distance=True)

		plain = np.atleast_2d(plain) if plain.size else np.empty((0, 2))
		assert len(with_distance) == plain.shape[0]
		assert np.array_equal(np.sort(with_distance["ra"]), np.sort(plain[:,0]))
