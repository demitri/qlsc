
import pytest
import numpy as np

from qlsc import QLSC, QLSCIndex

# QLSCIndex offers a dict-like surface over its keys: len(), "in", iteration over
# keys, and remove_point(). These tests cover both the in-memory index (whose
# connection must stay open for the life of the object) and a file-backed one
# (whose connections are opened and closed per operation).


def _populate(index):
	index.add_points(ra=np.array([10., 11., 12., 13.]),
					 dec=np.array([30., 30., 30., 30.]),
					 keys=["a", "b", "b", "c"])
	index.add_point(20., 40.) # deliberately keyless
	return index


@pytest.fixture(params=["memory", "file"])
def index(request, tmp_path):
	'''
	A populated index, once in memory and once backed by a file.
	'''
	if request.param == "memory":
		index = QLSCIndex(qlsc=QLSC(depth=30))
	else:
		index = QLSCIndex(qlsc=QLSC(depth=30), filepath=tmp_path / "mapping.qlscidx")
	yield _populate(index)
	index.close()


def test_len(index):
	'''
	Test that len() counts every point, keyed or not.
	'''
	assert len(index) == 5
	assert len(index) == index.number_of_points


def test_len_empty_index():
	'''
	Test that an empty index has length zero - and is therefore falsy, as an
	empty dict is.
	'''
	with QLSCIndex(qlsc=QLSC(depth=30)) as empty:
		assert len(empty) == 0
		assert not empty

	with QLSCIndex(qlsc=QLSC(depth=30)) as filled:
		filled.add_point(1., 2.)
		assert filled


def test_contains(index):
	'''
	Test key membership.
	'''
	assert "a" in index
	assert "b" in index
	assert "c" in index
	assert "z" not in index
	assert "" not in index


def test_contains_rejects_non_strings(index):
	'''
	Test that a non-string key raises rather than silently reporting False.

	The key column has SQLite TEXT affinity, so a bound integer would be compared
	against its text representation - a mismatch that looks like an absent key.
	'''
	for bad in (5, None, 3.4, b"a", ["a"]):
		with pytest.raises(TypeError):
			bad in index


def test_keys(index):
	'''
	Test iterating the keys.

	Keys are not unique, so a key shared by two points is yielded twice; the
	keyless point is not represented at all.
	'''
	assert sorted(index.keys()) == ["a", "b", "b", "c"]
	assert sorted(index) == ["a", "b", "b", "c"]

	# every key yielded must also be reported by __contains__
	for key in index.keys():
		assert key in index


def test_keys_is_a_generator(index):
	'''
	Test that keys() returns a generator rather than a materialized list.
	'''
	keys = index.keys()
	assert iter(keys) is keys
	assert next(keys) is not None


def test_keys_of_empty_index():
	'''
	Test that an index with no points, and one with no keyed points, both
	yield nothing.
	'''
	with QLSCIndex(qlsc=QLSC(depth=30)) as empty:
		assert list(empty.keys()) == []
		empty.add_point(1., 2.)
		assert list(empty.keys()) == []
		assert len(empty) == 1


def test_abandoned_key_generator_does_not_break_the_index(tmp_path):
	'''
	Test that a partly-consumed keys() generator does not leave the index unusable.

	File-backed indices open a connection per operation; the generator holds one
	for the length of its iteration and must release it when it is closed.
	'''
	index = QLSCIndex(qlsc=QLSC(depth=30), filepath=tmp_path / "abandoned.qlscidx")
	_populate(index)

	keys = index.keys()
	next(keys)
	keys.close() # abandon it mid-iteration

	assert len(index) == 5
	assert sorted(index.keys()) == ["a", "b", "b", "c"]
	index.close()


def test_remove_point_by_key(index):
	'''
	Test removing by key: every point carrying the key goes.
	'''
	assert index.remove_point(key="b") == 2
	assert len(index) == 3
	assert "b" not in index
	assert sorted(index.keys()) == ["a", "c"]

	assert index.remove_point(key="a") == 1
	assert len(index) == 2


def test_remove_point_by_coordinate(index):
	'''
	Test removing by coordinate, including the keyless point (which cannot be
	removed any other way).
	'''
	assert index.remove_point(ra=20., dec=40.) == 1
	assert len(index) == 4
	assert sorted(index.keys()) == ["a", "b", "b", "c"]

	assert index.remove_point(ra=11., dec=30.) == 1
	assert len(index) == 3
	assert sorted(index.keys()) == ["a", "b", "c"]


def test_removed_points_are_gone_from_queries(index):
	'''
	Test that a removed point is no longer returned by a radial query - the
	point of removing it.
	'''
	before = index.radial_query(ra=11., dec=30., radius=0.1)
	assert before.shape == (2,)

	index.remove_point(ra=11., dec=30.)

	after = index.radial_query(ra=11., dec=30., radius=0.1)
	assert after.size == 0


def test_remove_point_normalizes_coordinates():
	'''
	Test that removing by coordinate applies the same normalization add_point
	does, so a point added with an out-of-range dec can be removed with the
	same values it was added with.
	'''
	with QLSCIndex(qlsc=QLSC(depth=30)) as index:
		index.add_point(10., 100.) # normalized on the way in
		assert len(index) == 1
		assert index.remove_point(ra=10., dec=100.) == 1
		assert len(index) == 0


def test_remove_missing_point_raises(index):
	'''
	Test that removing something that is not there raises KeyError rather than
	quietly reporting success.
	'''
	with pytest.raises(KeyError):
		index.remove_point(key="not-a-key")
	with pytest.raises(KeyError):
		index.remove_point(ra=123.4, dec=-56.7)

	# a near-miss coordinate is still a miss: no tolerance is applied
	with pytest.raises(KeyError):
		index.remove_point(ra=10.0000001, dec=30.)

	assert len(index) == 5 # nothing was removed by the failed attempts


def test_remove_point_argument_validation(index):
	'''
	Test that ambiguous or incomplete arguments raise ValueError.
	'''
	with pytest.raises(ValueError):
		index.remove_point()
	with pytest.raises(ValueError):
		index.remove_point(ra=10., dec=30., key="a")
	with pytest.raises(ValueError):
		index.remove_point(ra=10.)
	with pytest.raises(ValueError):
		index.remove_point(dec=30.)

	with pytest.raises(TypeError):
		index.remove_point(key=5)

	assert len(index) == 5


def test_remove_point_accepts_scalar_like_coordinates(index):
	'''
	Test that the coordinate form accepts NumPy scalars and 0-d/size-1 arrays,
	as every other coordinate-taking API does.
	'''
	assert index.remove_point(ra=np.float64(10.), dec=np.array(30.)) == 1
	assert index.remove_point(ra=np.array([11.]), dec=np.array([30.])) == 1
	assert len(index) == 3


def test_removal_survives_reopening(tmp_path):
	'''
	Test that a removal is committed to a file-backed index, not just to the
	connection that performed it.
	'''
	path = tmp_path / "removal.qlscidx"

	index = QLSCIndex(qlsc=QLSC(depth=30), filepath=path)
	_populate(index)
	index.remove_point(key="b")
	index.close()

	reopened = QLSCIndex(qlsc=QLSC(depth=30), filepath=path)
	assert len(reopened) == 3
	assert "b" not in reopened
	assert sorted(reopened.keys()) == ["a", "c"]
	reopened.close()


def test_mapping_operations_on_a_closed_index():
	'''
	Test that every new operation refuses to work on a closed index rather than
	silently opening a fresh, empty database.
	'''
	index = _populate(QLSCIndex(qlsc=QLSC(depth=30)))
	index.close()

	with pytest.raises(RuntimeError):
		len(index)
	with pytest.raises(RuntimeError):
		"a" in index
	with pytest.raises(RuntimeError):
		index.keys()
	with pytest.raises(RuntimeError):
		list(index)
	with pytest.raises(RuntimeError):
		index.remove_point(key="a")
	with pytest.raises(RuntimeError):
		index.remove_point(ra=10., dec=30.)


def test_add_after_remove(index):
	'''
	Test that a removed coordinate can be added again (the table's UNIQUE(ra,dec)
	constraint and the INSERT OR IGNORE would otherwise make removal one-way).
	'''
	index.remove_point(ra=10., dec=30.)
	assert len(index) == 4

	index.add_point(10., 30., key="a-again")
	assert len(index) == 5
	assert "a-again" in index
