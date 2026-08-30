
import math

import pytest
import numpy as np

import qlsc
from qlsc import QLSC, QLSCIndex, crossmatch, sindist
from qlsc.generate import sunflower_points_on_sphere

# qlsc.crossmatch spatially joins two catalogs. Its correctness is checked against a
# brute-force O(n*m) sweep of every possible pair through the same exact predicate the
# index filters with (sindist), which shares no code with the pixelization: if the ipix
# range search ever under-covered a cone, these tests would see the missing pairs.


def _sindist_radius(radius_deg):
	''' The sin^2(d/2) threshold q3c's sindist is compared against. '''
	return math.pow(math.sin(math.radians(radius_deg)/2.), 2)


def _brute_force_pairs(points_a, points_b, radius_deg):
	'''
	Every pair within the radius, found by comparing all n*m combinations.

	Coordinates that are equivalent representations of the same point on the sphere
	(ra=360 vs ra=0, any ra at a pole) are treated as coincident, as the index does;
	sindist leaves a tiny nonzero residual between them which a zero or near-zero radius
	would otherwise reject.
	'''
	threshold = _sindist_radius(radius_deg)
	pairs = set()
	for i,(ra_a, dec_a) in enumerate(points_a):
		for j,(ra_b, dec_b) in enumerate(points_b):
			same_point = (abs(dec_a) == 90 and dec_a == dec_b) or (ra_a % 360. == ra_b % 360. and dec_a == dec_b)
			if same_point or (radius_deg > 0 and sindist(ra_a, dec_a, ra_b, dec_b) <= threshold):
				pairs.add((i, j))
	return pairs


def _index(points, depth=30, keys=None, filepath=None):
	index = QLSCIndex(qlsc=QLSC(depth=depth), filepath=filepath)
	points = np.asarray(points, dtype=np.double)
	index.add_points(ra=points[:,0], dec=points[:,1], keys=keys)
	return index


def _pair_indices(result, points_a, points_b):
	'''
	Map the returned coordinate pairs back to positions in the input point lists.
	'''
	def locate(points, ra, dec):
		for n,(point_ra, point_dec) in enumerate(points):
			if math.isclose(point_ra, ra, abs_tol=1e-12) and math.isclose(point_dec, dec, abs_tol=1e-12):
				return n
		raise AssertionError(f"the returned point (ra={ra}, dec={dec}) is not in the catalog it came from")

	return set((locate(points_a, row["ra_a"], row["dec_a"]), locate(points_b, row["ra_b"], row["dec_b"]))
			   for row in result)


def _rotated(points, d_ra, d_dec):
	'''
	Shift a set of points in ra and dec, wrapping ra and reflecting dec at the poles, to
	put a distribution somewhere else on the sky.
	'''
	moved = np.array(points, dtype=np.double, copy=True)
	moved[:,0] = (moved[:,0] + d_ra) % 360.
	moved[:,1] = np.clip(moved[:,1] + d_dec, -90., 90.)
	return moved


# ---------------------------------------------------------------------------
# agreement with the brute-force oracle
# ---------------------------------------------------------------------------

UNIFORM_A = sunflower_points_on_sphere(300)

# UNIFORM_B is what a real cross-match sees: an all-sky catalog of its own, of which the
# first 150 entries are counterparts of UNIFORM_A displaced by a few tenths of a degree,
# so that a modest radius finds a substantial number of pairs and a large one finds more.
_uniform_rng = np.random.default_rng(1414)
UNIFORM_B = _rotated(sunflower_points_on_sphere(400), d_ra=0.37, d_dec=0.11)
UNIFORM_B[:150] = _rotated(UNIFORM_A[:150], d_ra=0., d_dec=0.) + \
				  _uniform_rng.normal(scale=0.15, size=(150, 2))
UNIFORM_B[:,0] %= 360.
UNIFORM_B[:,1] = np.clip(UNIFORM_B[:,1], -90., 90.)

# a clustered pair of catalogs: many points inside a few degrees, so most cones hold
# several counterparts and nearest_only has something to choose between
_rng = np.random.default_rng(2026)
CLUSTERED_A = np.column_stack((_rng.normal(30., 2.0, 250) % 360., np.clip(_rng.normal(-15., 2.0, 250), -90, 90)))
CLUSTERED_B = np.column_stack((_rng.normal(30., 2.0, 300) % 360., np.clip(_rng.normal(-15., 2.0, 300), -90, 90)))

# the poles and the ra=0 seam, where equivalent representations of one point exist
SEAM_A = np.array([[0., 90.], [10., 89.99], [359.99, 0.], [0.01, 0.], [0., -90.], [180., -89.99]])
SEAM_B = np.array([[123., 90.], [10.005, 89.99], [0.005, 0.], [359.995, 0.], [77., -90.], [180.01, -89.985]])


@pytest.mark.parametrize("radius", [0.0, 0.01, 0.5, 3.0])
@pytest.mark.parametrize("name", ["uniform", "clustered", "seam"])
def test_agrees_with_brute_force(name, radius):
	'''
	Test that the pairs found are exactly the pairs a full n*m sweep finds.
	'''
	points_a, points_b = {"uniform"  : (UNIFORM_A, UNIFORM_B),
						  "clustered": (CLUSTERED_A, CLUSTERED_B),
						  "seam"     : (SEAM_A, SEAM_B)}[name]

	index_a = _index(points_a)
	index_b = _index(points_b)

	result = crossmatch(index_a, index_b, radius=radius)
	expected = _brute_force_pairs(points_a, points_b, radius)

	assert _pair_indices(result, points_a, points_b) == expected

	# guard against a vacuous comparison: two empty sets compare equal, so the cases that
	# are supposed to find something must be seen to find it
	if radius >= 0.5:
		assert len(expected) > 0, "this case is meant to exercise real matches"

	index_a.close()
	index_b.close()


def test_separations_are_correct():
	'''
	Test that the reported separation is the angular separation of the pair.
	'''
	index_a = _index(UNIFORM_A)
	index_b = _index(UNIFORM_B)

	result = crossmatch(index_a, index_b, radius=2.0)
	assert len(result) > 0

	for row in result:
		expected = math.degrees(2. * math.asin(math.sqrt(
			min(sindist(row["ra_a"], row["dec_a"], row["ra_b"], row["dec_b"]), 1.))))
		assert math.isclose(row["distance"], expected, abs_tol=1e-9)
		assert row["distance"] <= 2.0 + 1e-9

	index_a.close()
	index_b.close()


def test_iteration_direction_does_not_change_the_answer():
	'''
	Test that the smaller catalog being A or B gives the same pairs in the same order.

	The search iterates over whichever index is smaller; the output ordering is defined
	by the data, so swapping which side is bigger must not change the returned table.
	'''
	small = UNIFORM_A[:80]
	large = UNIFORM_B

	index_small = _index(small)
	index_large = _index(large)

	forward = crossmatch(index_small, index_large, radius=3.0)  # A is the smaller: iterates over A
	repeated = crossmatch(index_small, index_large, radius=3.0)
	swapped = crossmatch(index_large, index_small, radius=3.0)  # A is the larger: iterates over B

	assert len(forward) > 0
	assert np.array_equal(forward["ra_a"], repeated["ra_a"]), "the ordering must be stable across calls"
	assert len(forward) == len(swapped), "the relation is symmetric; only the roles are swapped"

	forward_pairs = set(zip(np.round(forward["ra_a"], 9), np.round(forward["ra_b"], 9)))
	swapped_pairs = set(zip(np.round(swapped["ra_b"], 9), np.round(swapped["ra_a"], 9)))
	assert forward_pairs == swapped_pairs

	index_small.close()
	index_large.close()


# ---------------------------------------------------------------------------
# nearest_only
# ---------------------------------------------------------------------------

def test_nearest_only_keeps_the_closest_counterpart():
	'''
	Test that nearest_only returns one row per index_a point, and that the row it keeps
	is the closest of the pairs the full match returns.
	'''
	index_a = _index(CLUSTERED_A)
	index_b = _index(CLUSTERED_B)

	everything = crossmatch(index_a, index_b, radius=0.5)
	nearest = crossmatch(index_a, index_b, radius=0.5, nearest_only=True)

	all_pairs = dict()
	for row in everything:
		all_pairs.setdefault((row["ra_a"], row["dec_a"]), list()).append((row["distance"], row["ra_b"], row["dec_b"]))

	assert len(nearest) == len(all_pairs), "one row per index_a point that matched anything"
	assert len(nearest) < len(everything), "this catalog must have points with several counterparts"

	seen = set()
	for row in nearest:
		point = (row["ra_a"], row["dec_a"])
		assert point not in seen, "a point of index_a may appear only once"
		seen.add(point)

		best = min(all_pairs[point])
		assert math.isclose(row["distance"], best[0], abs_tol=1e-12)
		assert (row["ra_b"], row["dec_b"]) == (best[1], best[2])

	index_a.close()
	index_b.close()


def test_nearest_only_does_not_invent_rows():
	'''
	Test that a point of index_a with no counterpart stays absent under nearest_only.
	'''
	points_a = [[10., 30.], [200., -70.]]
	points_b = [[10.001, 30.]]

	index_a = _index(points_a)
	index_b = _index(points_b)

	nearest = crossmatch(index_a, index_b, radius=0.1, nearest_only=True)

	assert len(nearest) == 1
	assert math.isclose(nearest[0]["ra_a"], 10.)

	index_a.close()
	index_b.close()


def test_nearest_only_tie_break_is_deterministic():
	'''
	Test that two counterparts at exactly the same separation resolve the same way every
	time, by ascending (ra_b, dec_b) as documented.
	'''
	index_a = _index([[10., 0.]])
	index_b = _index([[10.1, 0.], [9.9, 0.]])

	for _ in range(3):
		nearest = crossmatch(index_a, index_b, radius=0.5, nearest_only=True)
		assert len(nearest) == 1
		assert math.isclose(nearest[0]["ra_b"], 9.9), "the lower ra_b wins the tie"

	index_a.close()
	index_b.close()


# ---------------------------------------------------------------------------
# keys, ordering, shapes
# ---------------------------------------------------------------------------

def test_keys_are_returned():
	'''
	Test that the keys of both catalogs come back with the pairs.
	'''
	keys_a = [f"a{n}" for n in range(len(UNIFORM_A))]
	keys_b = [f"b{n}" for n in range(len(UNIFORM_B))]

	index_a = _index(UNIFORM_A, keys=keys_a)
	index_b = _index(UNIFORM_B, keys=keys_b)

	result = crossmatch(index_a, index_b, radius=2.0)

	assert result.dtype.names == ("ra_a", "dec_a", "key_a", "ra_b", "dec_b", "key_b", "distance")
	assert len(result) > 0
	assert all(key.startswith("a") for key in result["key_a"])
	assert all(key.startswith("b") for key in result["key_b"])

	# the key must belong to the point it is reported with
	lookup_a = {(round(ra, 9), round(dec, 9)): key for (ra,dec),key in zip(UNIFORM_A, keys_a)}
	for row in result:
		assert lookup_a[(round(row["ra_a"], 9), round(row["dec_a"], 9))] == row["key_a"]

	index_a.close()
	index_b.close()


def test_points_without_keys_report_none():
	'''
	Test that a point stored without a key reports None rather than the string "None".
	'''
	index_a = _index([[10., 30.]], keys=["only-a"])
	index_b = _index([[10.001, 30.]])

	result = crossmatch(index_a, index_b, radius=0.1)

	assert len(result) == 1
	assert result[0]["key_a"] == "only-a"
	assert result[0]["key_b"] is None

	index_a.close()
	index_b.close()


def test_return_keys_false_omits_the_fields():
	'''
	Test that the key fields can be left out.
	'''
	index_a = _index(UNIFORM_A[:50])
	index_b = _index(UNIFORM_B[:50])

	result = crossmatch(index_a, index_b, radius=5.0, return_keys=False)

	assert result.dtype.names == ("ra_a", "dec_a", "ra_b", "dec_b", "distance")

	index_a.close()
	index_b.close()


def test_ordering_is_grouped_and_sorted():
	'''
	Test the documented ordering: grouped by (ra_a, dec_a) ascending, then by increasing
	separation within a group.
	'''
	index_a = _index(CLUSTERED_A)
	index_b = _index(CLUSTERED_B)

	result = crossmatch(index_a, index_b, radius=0.6)
	assert len(result) > 0

	previous_point = None
	previous_distance = None
	seen_points = set()
	for row in result:
		point = (row["ra_a"], row["dec_a"])
		if point != previous_point:
			assert point not in seen_points, "the rows of one index_a point must be contiguous"
			if previous_point is not None:
				assert point > previous_point, "groups must be in ascending (ra_a, dec_a) order"
			seen_points.add(point)
			previous_point = point
			previous_distance = row["distance"]
		else:
			assert row["distance"] >= previous_distance - 1e-12, "a group must be sorted by separation"
			previous_distance = row["distance"]

	index_a.close()
	index_b.close()


# ---------------------------------------------------------------------------
# empty results, self-matches, mixed depths, file-backed indices
# ---------------------------------------------------------------------------

def test_no_pairs_within_the_radius():
	'''
	Test that two catalogs on opposite sides of the sky return an empty, correctly typed
	table rather than failing.
	'''
	index_a = _index([[10., 30.], [11., 30.]])
	index_b = _index([[190., -30.], [191., -30.]])

	result = crossmatch(index_a, index_b, radius=1.0)

	assert len(result) == 0
	assert result.dtype.names == ("ra_a", "dec_a", "key_a", "ra_b", "dec_b", "key_b", "distance")

	index_a.close()
	index_b.close()


def test_empty_index():
	'''
	Test that an empty catalog on either side gives an empty result.
	'''
	populated = _index(UNIFORM_A[:20])
	empty = QLSCIndex(qlsc=QLSC(depth=30))

	assert len(crossmatch(populated, empty, radius=10.0)) == 0
	assert len(crossmatch(empty, populated, radius=10.0)) == 0
	assert len(crossmatch(empty, empty, radius=10.0)) == 0

	populated.close()
	empty.close()


def test_self_match():
	'''
	Test cross-matching an index with itself: every point pairs with itself at zero
	separation, plus any neighbours within the radius.
	'''
	index = _index(CLUSTERED_A)

	result = crossmatch(index, index, radius=0.3)
	assert _pair_indices(result, CLUSTERED_A, CLUSTERED_A) == \
		   _brute_force_pairs(CLUSTERED_A, CLUSTERED_A, 0.3)

	# every point is paired with itself
	self_pairs = [row for row in result if row["distance"] == 0.0]
	assert len(self_pairs) == len(CLUSTERED_A)
	for row in self_pairs:
		assert (row["ra_a"], row["dec_a"]) == (row["ra_b"], row["dec_b"])

	# and under nearest_only the self-pair is what survives, as documented
	nearest = crossmatch(index, index, radius=0.3, nearest_only=True)
	assert len(nearest) == len(CLUSTERED_A)
	assert (nearest["distance"] == 0.0).all()

	index.close()


def test_self_match_at_zero_radius():
	'''
	Test a zero-radius self-match: exactly the self-pairs, nothing else.
	'''
	index = _index(UNIFORM_A[:60])

	result = crossmatch(index, index, radius=0.0)

	assert len(result) == 60
	assert (result["distance"] == 0.0).all()
	assert np.array_equal(result["ra_a"], result["ra_b"])

	index.close()


@pytest.mark.parametrize("depth_a, depth_b", [(30, 30), (30, 8), (8, 30), (5, 12), (0, 30), (30, 0)])
def test_mixed_depths(depth_a, depth_b):
	'''
	Test that the two indices do not have to share a QLSC depth: each search uses the
	depth of the index being searched, and the coarse ranges are a superset either way.
	'''
	points_a, points_b = UNIFORM_A[:120], UNIFORM_B[:150]

	index_a = _index(points_a, depth=depth_a)
	index_b = _index(points_b, depth=depth_b)

	result = crossmatch(index_a, index_b, radius=2.0)

	assert _pair_indices(result, points_a, points_b) == _brute_force_pairs(points_a, points_b, 2.0)

	index_a.close()
	index_b.close()


def test_file_backed_indices(tmp_path):
	'''
	Test cross-matching two indices stored on disk, and one of each kind.
	'''
	points_a, points_b = UNIFORM_A[:100], UNIFORM_B[:120]
	expected = _brute_force_pairs(points_a, points_b, 3.0)

	file_a = _index(points_a, filepath=tmp_path / "a.qlscidx")
	file_b = _index(points_b, filepath=tmp_path / "b.qlscidx")
	memory_b = _index(points_b)

	assert _pair_indices(crossmatch(file_a, file_b, radius=3.0), points_a, points_b) == expected
	assert _pair_indices(crossmatch(file_a, memory_b, radius=3.0), points_a, points_b) == expected

	file_a.close()
	file_b.close()
	memory_b.close()


def test_astropy_quantity_radius():
	'''
	Test that the radius may be an astropy Quantity.
	'''
	astropy = pytest.importorskip("astropy.units")

	index_a = _index(UNIFORM_A[:60])
	index_b = _index(UNIFORM_B[:80])

	as_quantity = crossmatch(index_a, index_b, radius=120 * astropy.arcmin)
	as_degrees = crossmatch(index_a, index_b, radius=2.0)

	assert np.array_equal(as_quantity["ra_a"], as_degrees["ra_a"])
	assert np.array_equal(as_quantity["ra_b"], as_degrees["ra_b"])

	index_a.close()
	index_b.close()


# ---------------------------------------------------------------------------
# invalid input
# ---------------------------------------------------------------------------

def test_missing_arguments_are_refused():
	index = _index([[10., 30.]])

	with pytest.raises(ValueError, match="index_a"):
		crossmatch(None, index, radius=1.0)
	with pytest.raises(ValueError, match="index_b"):
		crossmatch(index, None, radius=1.0)
	with pytest.raises(ValueError, match="radius"):
		crossmatch(index, index)

	index.close()


def test_non_index_argument_is_refused():
	index = _index([[10., 30.]])

	with pytest.raises(TypeError, match="QLSCIndex"):
		crossmatch(index, QLSC(depth=30), radius=1.0)

	index.close()


def test_out_of_range_radius_is_refused():
	index = _index([[10., 30.]])

	with pytest.raises(ValueError, match=r"\[0,180\]"):
		crossmatch(index, index, radius=181.)
	with pytest.raises(ValueError, match=r"\[0,180\]"):
		crossmatch(index, index, radius=-1.)
	with pytest.raises(ValueError, match="finite"):
		crossmatch(index, index, radius=float("nan"))

	index.close()


def test_closed_index_is_refused():
	index_a = _index([[10., 30.]])
	index_b = _index([[10., 30.]])
	index_b.close()

	with pytest.raises(RuntimeError, match="closed"):
		crossmatch(index_a, index_b, radius=1.0)

	index_a.close()


def test_crossmatch_is_exported():
	'''
	Test that the function is reachable as qlsc.crossmatch, which is how it is documented.
	'''
	assert qlsc.crossmatch is crossmatch


def test_all_points_reads_the_whole_index():
	'''
	Test the private read behind crossmatch, including the paths taken when the row count
	it sizes its arrays from turns out to be wrong.

	The count and the read are two statements, and this class opens a connection per
	operation precisely so that it can be used from several threads, so the index can
	change between them. No point may be dropped to make the answer fit the count.
	'''
	from unittest import mock

	points = UNIFORM_A[:40]
	keys = [f"k{n}" for n in range(len(points))]
	index = _index(points, keys=keys)

	ra, dec, read_keys = index._all_points(with_keys=True)
	assert len(ra) == len(points)
	assert set(read_keys) == set(keys)
	assert set(np.round(ra, 9)) == set(np.round(points[:,0], 9))

	# without keys
	ra_only, dec_only, no_keys = index._all_points()
	assert no_keys is None
	assert np.array_equal(np.sort(ra_only), np.sort(ra))

	# the count under-reports (the index "grew"): every point must still be read
	with mock.patch.object(type(index), "number_of_points", new_callable=mock.PropertyMock) as count:
		count.return_value = 3
		grown_ra, grown_dec, grown_keys = index._all_points(with_keys=True)
	assert len(grown_ra) == len(points) and len(grown_keys) == len(points)
	assert np.array_equal(np.sort(grown_ra), np.sort(ra))
	assert set(grown_keys) == set(keys)

	# the count over-reports (the index "shrank"): the arrays must be trimmed, not padded
	with mock.patch.object(type(index), "number_of_points", new_callable=mock.PropertyMock) as count:
		count.return_value = 500
		trimmed_ra, trimmed_dec, trimmed_keys = index._all_points(with_keys=True)
	assert len(trimmed_ra) == len(points) and len(trimmed_keys) == len(points)
	assert np.array_equal(np.sort(trimmed_ra), np.sort(ra))

	index.close()


def test_all_points_of_an_empty_index():
	'''
	Test that reading an empty index gives empty arrays of the right dtype.
	'''
	index = QLSCIndex(qlsc=QLSC(depth=30))
	ra, dec, keys = index._all_points(with_keys=True)

	assert ra.shape == (0,) and ra.dtype == np.double
	assert dec.shape == (0,) and dec.dtype == np.double
	assert keys == []

	index.close()
