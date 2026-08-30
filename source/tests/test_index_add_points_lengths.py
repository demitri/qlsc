
import pytest
import numpy as np

from qlsc import QLSC, QLSCIndex
from qlsc.qlsc_base import _zip_exact

# add_points builds its rows with zip(), which stops at its shortest input. Every
# length mismatch must therefore be an error: a short 'keys' list previously stored
# only as many points as there were keys and dropped the rest without a word.


def _index():
	return QLSCIndex(qlsc=QLSC(depth=30))


@pytest.mark.parametrize("n_coordinates, n_keys", [(3, 2), (2, 3), (5, 1), (1, 5), (4, 0)])
def test_mismatched_key_count_raises(n_coordinates, n_keys):
	'''
	Test that a key list of the wrong length raises instead of truncating the insert.

	Regression test: the check was written as the chained comparison
	"len(ra) != len(dec) != len(keys)", which Python evaluates as
	(len(ra) != len(dec)) and (len(dec) != len(keys)) - False whenever the
	coordinates agree, whatever the key count is. Three coordinates and two keys
	stored two points.
	'''
	ra = np.linspace(10., 20., n_coordinates)
	dec = np.linspace(-5., 5., n_coordinates)
	keys = [f"k{i}" for i in range(n_keys)]

	with _index() as index:
		with pytest.raises(ValueError):
			index.add_points(ra=ra, dec=dec, keys=keys)
		assert len(index) == 0, "a rejected call must not have stored anything"


@pytest.mark.parametrize("n_points, n_keys", [(3, 2), (2, 3)])
def test_mismatched_key_count_raises_for_the_points_form(n_points, n_keys):
	'''
	Test the same for the 'points' argument form (whose check was already correct).
	'''
	points = np.stack([np.linspace(10., 20., n_points), np.linspace(-5., 5., n_points)], axis=1)
	keys = [f"k{i}" for i in range(n_keys)]

	with _index() as index:
		with pytest.raises(ValueError):
			index.add_points(points=points, keys=keys)
		assert len(index) == 0


def test_matching_lengths_still_work():
	'''
	Test that the guard does not reject correct input.
	'''
	ra = np.array([10., 11., 12.])
	dec = np.array([30., 31., 32.])

	with _index() as index:
		index.add_points(ra=ra, dec=dec, keys=["a", "b", "c"])
		assert len(index) == 3
		assert sorted(index.keys()) == ["a", "b", "c"]

	with _index() as index:
		index.add_points(points=np.stack([ra, dec], axis=1), keys=["a", "b", "c"])
		assert len(index) == 3


def test_mismatched_coordinate_lengths_raise():
	'''
	Test that unequal ra and dec arrays raise, with and without keys.
	'''
	with _index() as index:
		with pytest.raises(ValueError):
			index.add_points(ra=np.array([1., 2., 3.]), dec=np.array([1., 2.]))
		with pytest.raises(ValueError):
			index.add_points(ra=np.array([1., 2., 3.]), dec=np.array([1., 2.]), keys=["a", "b", "c"])
		assert len(index) == 0


def test_both_argument_forms_at_once_raises():
	'''
	Test that giving both 'points' and 'ra','dec' raises.

	Regression test: the check was "all([x is True for x in [ra, dec, points]])",
	which no array ever satisfies, so 'points' silently won and the ra,dec pair was
	discarded.
	'''
	points = np.array([[10., 30.], [11., 31.]])

	with _index() as index:
		with pytest.raises(ValueError):
			index.add_points(ra=np.array([50., 51.]), dec=np.array([-20., -21.]), points=points)
		assert len(index) == 0


def test_zip_exact_guard():
	'''
	Test the guard itself, which is what enforces the invariant at the point the
	rows are built rather than only where the arguments are validated.
	'''
	assert list(_zip_exact([1, 2], "ab")) == [(1, "a"), (2, "b")]
	assert list(_zip_exact([], [])) == []

	with pytest.raises(ValueError):
		_zip_exact([1, 2, 3], [1, 2])
	with pytest.raises(ValueError):
		_zip_exact([1, 2], [1, 2], [1])
