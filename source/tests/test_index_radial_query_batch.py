
import math

import pytest
import numpy as np

from qlsc import QLSC, QLSCIndex
from qlsc.generate import sunflower_points_on_sphere

# QLSCIndex.radial_query accepts a list of targets as well as a single coordinate, for
# matching a target list against the index. The contract these tests pin down:
#
#   * a batch query is one whose ra/dec have a shape AND a size other than 1, or one
#     with return_counts set; a size-1 array is a scalar, as it is everywhere in qlsc
#   * the result is always a recarray whose first field is "target_index"
#   * a target that matches nothing contributes no rows; return_counts makes it explicit
#   * the answer is the same as calling radial_query once per target

CATALOG = sunflower_points_on_sphere(4000)


def _index(points=CATALOG, depth=30, keys=None):
	index = QLSCIndex(qlsc=QLSC(depth=depth))
	points = np.asarray(points, dtype=np.double)
	index.add_points(ra=points[:,0], dec=points[:,1], keys=keys)
	return index


def _targets(n, seed, catalog=CATALOG, jitter=0.0):
	'''
	Pick n catalog positions at random, optionally moved off the catalog point.
	'''
	rng = np.random.default_rng(seed)
	chosen = catalog[rng.choice(len(catalog), n, replace=False)]
	if jitter:
		chosen = chosen + rng.normal(scale=jitter, size=chosen.shape)
	return chosen[:,0], chosen[:,1]


def _as_pair_set(matches):
	return set(zip(np.round(matches["ra"], 9), np.round(matches["dec"], 9)))


@pytest.mark.parametrize("radius", [0.5, 2.0, 10.0])
@pytest.mark.parametrize("n_targets", [1, 7, 40])
def test_batch_agrees_with_individual_queries(radius, n_targets):
	'''
	Test that a batch query returns exactly what the same targets return one at a time.
	'''
	index = _index()
	ra, dec = _targets(n_targets, seed=11 + n_targets, jitter=0.3)

	matches, counts = index.radial_query(ra=ra, dec=dec, radius=radius,
										 return_counts=True, return_distance=True)

	assert counts.dtype == np.int64
	assert counts.shape == (n_targets,)
	assert counts.sum() == len(matches)

	for n in range(n_targets):
		expected = index.radial_query(float(ra[n]), float(dec[n]), radius, return_distance=True)
		mine = matches[matches["target_index"] == n]

		assert len(mine) == counts[n]
		assert _as_pair_set(mine) == _as_pair_set(expected)
		assert np.allclose(np.sort(mine["distance"]), np.sort(expected["distance"]))

	index.close()


def test_matches_are_grouped_in_target_order():
	'''
	Test that the flat table is ordered by target, in the order the targets were given.
	'''
	index = _index()
	ra, dec = _targets(12, seed=5, jitter=0.5)

	matches = index.radial_query(ra=ra, dec=dec, radius=3.0)

	assert (np.diff(matches["target_index"]) >= 0).all(), "matches must be grouped by target"
	assert matches["target_index"].min() >= 0
	assert matches["target_index"].max() < 12

	index.close()


def test_zero_match_targets_are_counted_not_dropped():
	'''
	Test the documented representation of a target that matches nothing: no rows in the
	flat table, and a zero in the counts array.

	The middle target here is deliberately placed in empty sky with a tiny radius.
	'''
	index = _index()
	catalog_ra, catalog_dec = CATALOG[10,0], CATALOG[10,1]

	ra = np.array([catalog_ra, catalog_ra + 0.4, catalog_ra])
	dec = np.array([catalog_dec, catalog_dec, catalog_dec])

	matches, counts = index.radial_query(ra=ra, dec=dec, radius=1e-7, return_counts=True)

	assert counts[0] == 1 and counts[2] == 1
	assert counts[1] == 0, "the target in empty sky must be represented by a zero count"
	assert 1 not in set(matches["target_index"]), "a zero-match target contributes no rows"
	assert len(matches) == counts.sum()

	# the same information, reconstructed the way the docstring describes
	assert np.array_equal(np.bincount(matches["target_index"], minlength=len(ra)), counts)

	index.close()


def test_all_targets_matching_nothing():
	'''
	Test a batch in which every target matches nothing: an empty table, all-zero counts.
	'''
	index = _index(points=[[10., 30.], [11., 30.]])

	matches, counts = index.radial_query(ra=[200., 210.], dec=[-40., -50.], radius=1.0,
										 return_counts=True)

	assert len(matches) == 0
	assert matches.dtype.names[0] == "target_index"
	assert matches.dtype["target_index"] == np.int64, \
		"an empty result must keep an integer target_index, so it can still subscript an array"
	assert np.array_equal(counts, [0, 0])

	index.close()


def test_empty_target_list():
	'''
	Test that a target list with no targets is an empty batch, not an error and not a
	scalar query.
	'''
	index = _index(points=[[10., 30.]])

	matches, counts = index.radial_query(ra=[], dec=[], radius=1.0, return_counts=True)

	assert len(matches) == 0
	assert counts.shape == (0,)
	assert "target_index" in matches.dtype.names

	index.close()


def test_per_target_radius():
	'''
	Test that a radius array applies one radius to each target.
	'''
	index = _index()
	ra, dec = _targets(6, seed=3, jitter=0.2)
	radii = np.array([0.2, 1.0, 2.0, 4.0, 8.0, 0.5])

	matches, counts = index.radial_query(ra=ra, dec=dec, radius=radii, return_counts=True)

	for n in range(len(ra)):
		expected = index.radial_query(float(ra[n]), float(dec[n]), float(radii[n]))
		assert counts[n] == (0 if expected.size == 0 else max(expected.size // 2, 1))
		assert _as_pair_set(matches[matches["target_index"] == n]) == \
			   set(zip(np.round(np.atleast_2d(expected)[:,0], 9), np.round(np.atleast_2d(expected)[:,1], 9)))

	index.close()


def test_keys_and_distances_together():
	'''
	Test the full set of fields, and their order.
	'''
	keys = [f"src-{n}" for n in range(len(CATALOG))]
	index = _index(keys=keys)
	ra, dec = _targets(5, seed=19, jitter=0.4)

	matches = index.radial_query(ra=ra, dec=dec, radius=2.0, return_key=True, return_distance=True)

	assert matches.dtype.names == ("target_index", "ra", "dec", "key", "distance")
	assert all(key.startswith("src-") for key in matches["key"])
	assert (matches["distance"] <= 2.0 + 1e-9).all()

	index.close()


def test_poles_and_seam_targets():
	'''
	Test targets at the poles and either side of ra=0, where equivalent coordinate
	representations of the same point have to match.
	'''
	points = [[0., 90.], [123., 90.], [0., -90.], [359.99, 0.], [0.01, 0.], [180., 45.]]
	index = _index(points=points)

	ra = np.array([200., 0., 0., 359.995, 180.])
	dec = np.array([90., -90., 0., 0., 45.])

	matches, counts = index.radial_query(ra=ra, dec=dec, radius=0.05, return_counts=True)

	for n in range(len(ra)):
		expected = index.radial_query(float(ra[n]), float(dec[n]), 0.05)
		assert counts[n] == (0 if expected.size == 0 else max(expected.size // 2, 1))

	# both stored representations of the north pole are found from any ra
	assert counts[0] == 2

	index.close()


def test_scalar_behaviour_is_unchanged():
	'''
	Test that scalar and scalar-like input still take the scalar path, unchanged.
	'''
	index = _index(points=[[10., 30.], [10.4, 30.], [200., -10.]])

	plain = index.radial_query(10., 30., 1.0)
	assert plain.dtype.names is None and plain.shape == (2, 2)

	# a size-1 array is a scalar value everywhere in qlsc, including here
	size_one = index.radial_query(ra=np.array([10.]), dec=np.array([30.]), radius=np.array([1.]))
	assert size_one.dtype.names is None and size_one.shape == (2, 2)

	zero_d = index.radial_query(ra=np.array(10.), dec=np.array(30.), radius=1.0)
	assert zero_d.dtype.names is None and zero_d.shape == (2, 2)

	index.close()


def test_single_target_list_via_return_counts():
	'''
	Test the documented way to ask for the batch contract with one target, which is
	otherwise indistinguishable from a scalar query.
	'''
	index = _index(points=[[10., 30.], [10.4, 30.]])

	matches, counts = index.radial_query(ra=[10.], dec=[30.], radius=1.0, return_counts=True)

	assert matches.dtype.names[0] == "target_index"
	assert set(matches["target_index"]) == {0}
	assert counts.shape == (1,) and counts[0] == 2

	index.close()


def test_file_backed_index(tmp_path):
	'''
	Test a batch query against a file-backed index, where a connection is opened and
	closed around every target's query.
	'''
	index = QLSCIndex(qlsc=QLSC(depth=20), filepath=tmp_path / "batch.qlscidx")
	index.add_points(ra=CATALOG[:500,0], dec=CATALOG[:500,1])

	ra, dec = CATALOG[:5,0], CATALOG[:5,1]
	matches, counts = index.radial_query(ra=ra, dec=dec, radius=5.0, return_counts=True)

	for n in range(5):
		assert counts[n] == len(index.radial_query(float(ra[n]), float(dec[n]), 5.0, return_distance=True))

	index.close()


def test_mismatched_target_lengths_are_refused():
	'''
	Test that ra and dec lists of different lengths raise.
	'''
	index = _index(points=[[10., 30.]])
	with pytest.raises(ValueError, match="does not match"):
		index.radial_query(ra=[1., 2., 3.], dec=[1., 2.], radius=1.0)
	index.close()


def test_mismatched_radius_length_is_refused():
	'''
	Test that a radius array of the wrong length raises rather than being broadcast.
	'''
	index = _index(points=[[10., 30.]])
	with pytest.raises(ValueError, match="does not match the number of targets"):
		index.radial_query(ra=[1., 2., 3.], dec=[1., 2., 3.], radius=[1., 2.])
	index.close()


def test_out_of_range_target_is_refused_before_any_query():
	'''
	Test that an invalid dec or radius anywhere in the list raises, naming the value.
	'''
	index = _index(points=[[10., 30.]])

	with pytest.raises(ValueError, match=r"\[-90,90\]"):
		index.radial_query(ra=[10., 20.], dec=[30., 91.], radius=1.0)

	with pytest.raises(ValueError, match=r"\[0,180\]"):
		index.radial_query(ra=[10., 20.], dec=[30., 40.], radius=[1., 200.])

	index.close()


def test_non_finite_target_is_refused():
	'''
	Test that a NaN target raises rather than reaching the C layer's undefined
	float-to-int conversions.
	'''
	index = _index(points=[[10., 30.]])
	with pytest.raises(ValueError, match="finite"):
		index.radial_query(ra=[10., np.nan], dec=[30., 40.], radius=1.0)
	index.close()


def test_two_dimensional_targets_are_refused():
	'''
	Test that a 2-D target array raises rather than being flattened.
	'''
	index = _index(points=[[10., 30.]])
	with pytest.raises(ValueError, match="one-dimensional"):
		index.radial_query(ra=np.zeros((2, 2)), dec=np.zeros((2, 2)), radius=1.0)
	index.close()


def test_astropy_quantity_radius():
	'''
	Test that an astropy Quantity radius works for a batch query too.
	'''
	astropy = pytest.importorskip("astropy.units")

	index = _index()
	ra, dec = _targets(4, seed=31, jitter=0.2)

	as_quantity, _ = index.radial_query(ra=ra, dec=dec, radius=60 * astropy.arcmin, return_counts=True)
	as_degrees, _ = index.radial_query(ra=ra, dec=dec, radius=1.0, return_counts=True)

	assert _as_pair_set(as_quantity) == _as_pair_set(as_degrees)

	index.close()


def test_astropy_quantity_coordinates():
	'''
	Test that Quantity coordinates are converted, not stripped of their units.

	np.asarray on a Quantity yields its numbers with the units discarded, so a target
	list given in radians would otherwise be read as degrees: the query would run at
	entirely the wrong place on the sky and report no matches rather than an error.
	'''
	astropy = pytest.importorskip("astropy.units")

	index = _index()
	ra, dec = _targets(6, seed=77, jitter=0.2)

	plain, plain_counts = index.radial_query(ra=ra, dec=dec, radius=1.5, return_counts=True)
	assert plain_counts.sum() > 0, "this case must find something for the comparison to mean anything"

	in_degrees, _ = index.radial_query(ra=ra * astropy.deg, dec=dec * astropy.deg,
									   radius=1.5, return_counts=True)
	in_radians, _ = index.radial_query(ra=np.radians(ra) * astropy.rad, dec=np.radians(dec) * astropy.rad,
									   radius=1.5, return_counts=True)

	assert _as_pair_set(in_degrees) == _as_pair_set(plain)
	assert _as_pair_set(in_radians) == _as_pair_set(plain), \
		"a radians Quantity must be converted, not read as degrees"

	index.close()


def test_astropy_quantity_scalar_coordinates():
	'''
	Test the same for the scalar path, which shares the conversion.
	'''
	astropy = pytest.importorskip("astropy.units")

	index = _index()
	ra, dec = float(CATALOG[100,0]), float(CATALOG[100,1])

	plain = index.radial_query(ra, dec, 2.0)
	assert plain.size > 0

	in_radians = index.radial_query(ra=math.radians(ra) * astropy.rad,
									dec=math.radians(dec) * astropy.rad, radius=2.0)

	assert np.allclose(np.atleast_2d(in_radians), np.atleast_2d(plain))

	index.close()
