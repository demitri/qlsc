
import pytest
import numpy as np

from qlsc import QLSC
from qlsc import q3c

# QLSC.ipix2ang, QLSC.ipix2ang_center and QLSC.ipix2xy accept either a single
# ipix value or an array-like of them, following the same scalar/array convention
# as ang2ipix: a scalar (or a 0-d array) returns scalars, anything with a shape
# returns arrays of that shape.

DEPTHS = [0, 1, 2, 5, 26, 30]


def _sample_ipix(q, count=25):
	'''
	Return a spread of valid ipix values for the given scheme.
	'''
	if q.nbins <= count:
		return np.arange(q.nbins, dtype=np.int64)
	# integer arithmetic only: np.linspace goes through float64, which rounds the
	# last value past nbins-1 at depth 30 (nbins = 6*2^60)
	step = (q.nbins - 1) // (count - 1)
	values = [i * step for i in range(count - 1)] + [q.nbins - 1]
	return np.array(values, dtype=np.int64)


@pytest.mark.parametrize("depth", DEPTHS)
def test_ipix2ang_array_matches_scalar(depth):
	'''
	Test that array-form ipix2ang returns exactly what the scalar form returns.
	'''
	q = QLSC(depth=depth)
	ipix = _sample_ipix(q)

	ra, dec = q.ipix2ang(ipix)

	assert isinstance(ra, np.ndarray) and isinstance(dec, np.ndarray)
	assert ra.shape == dec.shape == ipix.shape
	assert ra.dtype == np.double and dec.dtype == np.double

	for i, value in enumerate(ipix):
		scalar_ra, scalar_dec = q.ipix2ang(int(value))
		# the array path calls the identical C function; require bit equality
		assert ra[i] == scalar_ra
		assert dec[i] == scalar_dec


@pytest.mark.parametrize("depth", DEPTHS)
def test_ipix2xy_array_matches_scalar(depth):
	'''
	Test that array-form ipix2xy returns exactly what the scalar form returns.
	'''
	q = QLSC(depth=depth)
	ipix = _sample_ipix(q)

	facenum, x, y = q.ipix2xy(ipix)

	assert facenum.shape == x.shape == y.shape == ipix.shape
	assert np.issubdtype(facenum.dtype, np.integer)
	assert x.dtype == np.double and y.dtype == np.double

	for i, value in enumerate(ipix):
		scalar_facenum, scalar_x, scalar_y = q.ipix2xy(int(value))
		assert facenum[i] == scalar_facenum
		assert x[i] == scalar_x
		assert y[i] == scalar_y


@pytest.mark.parametrize("depth", DEPTHS)
def test_ipix2ang_center_array(depth):
	'''
	Test that the ipix2ang_center synonym vectorizes identically to ipix2ang.
	'''
	q = QLSC(depth=depth)
	ipix = _sample_ipix(q)

	ra, dec = q.ipix2ang(ipix)
	center_ra, center_dec = q.ipix2ang_center(ipix)

	assert np.array_equal(ra, center_ra)
	assert np.array_equal(dec, center_dec)


@pytest.mark.parametrize("depth", [2, 5, 26, 30])
def test_ipix2ang_array_round_trip(depth):
	'''
	Test that the pixel centers returned by array-form ipix2ang map back to the
	same ipix values through the (already vectorized) ang2ipix.
	'''
	q = QLSC(depth=depth)
	ipix = _sample_ipix(q, count=200)

	ra, dec = q.ipix2ang(ipix)
	assert np.array_equal(q.ang2ipix(ra=ra, dec=dec), ipix)


@pytest.mark.parametrize("method", ["ipix2ang", "ipix2ang_center", "ipix2xy"])
def test_scalar_like_input_returns_scalars(method):
	'''
	Test that scalar and 0-d array input returns scalars, not arrays.
	'''
	q = QLSC(depth=5)
	call = getattr(q, method)

	for value in (10, np.int64(10), np.int32(10), np.array(10)):
		result = call(value)
		for item in result:
			assert np.isscalar(item), f"{method}({value!r}) returned a non-scalar element: {item!r}"
		assert result == call(10)


@pytest.mark.parametrize("method", ["ipix2ang", "ipix2ang_center", "ipix2xy"])
def test_size_one_array_returns_arrays(method):
	'''
	Test that a length-1 array returns length-1 arrays rather than scalars.

	This is the same wart that ang2ipix had: a size-1 array satisfies the C
	scalar signature through a NumPy conversion that has since been deprecated.
	'''
	q = QLSC(depth=5)
	result = getattr(q, method)(np.array([10]))

	for item in result:
		assert isinstance(item, np.ndarray)
		assert item.shape == (1,)


@pytest.mark.parametrize("method", ["ipix2ang", "ipix2ang_center", "ipix2xy"])
def test_list_and_tuple_input(method):
	'''
	Test that plain Python sequences of integers are accepted.
	'''
	q = QLSC(depth=5)
	expected = getattr(q, method)(np.array([3, 17, 4095]))

	for form in ([3, 17, 4095], (3, 17, 4095)):
		result = getattr(q, method)(form)
		for a, b in zip(result, expected):
			assert np.array_equal(a, b)


@pytest.mark.parametrize("method", ["ipix2ang", "ipix2ang_center", "ipix2xy"])
def test_shape_is_preserved(method):
	'''
	Test that the shape of the input array is carried through to the outputs.
	'''
	q = QLSC(depth=5)
	ipix = np.arange(12, dtype=np.int64).reshape(3, 4)

	result = getattr(q, method)(ipix)
	for item in result:
		assert item.shape == (3, 4)

	flat = getattr(q, method)(ipix.ravel())
	for item, flat_item in zip(result, flat):
		assert np.array_equal(item.ravel(), flat_item)


@pytest.mark.parametrize("method", ["ipix2ang", "ipix2ang_center", "ipix2xy"])
def test_non_contiguous_input(method):
	'''
	Test that a strided (non-contiguous) array is handled correctly rather than
	being read as if it were contiguous.
	'''
	q = QLSC(depth=5)
	block = np.arange(0, 20, dtype=np.int64).reshape(10, 2)
	strided = block[:,1] # every other element

	result = getattr(q, method)(strided)
	expected = getattr(q, method)(np.ascontiguousarray(strided))

	for a, b in zip(result, expected):
		assert np.array_equal(a, b)


@pytest.mark.parametrize("method", ["ipix2ang", "ipix2ang_center", "ipix2xy"])
def test_integer_dtypes_accepted(method):
	'''
	Test that integer arrays of any width are accepted (ang2ipix returns int64,
	but catalog columns are often narrower).
	'''
	q = QLSC(depth=5)
	values = [3, 17, 4095]
	expected = getattr(q, method)(np.array(values, dtype=np.int64))

	for dtype in (np.int16, np.int32, np.uint16, np.uint32, np.uint64):
		result = getattr(q, method)(np.array(values, dtype=dtype))
		for a, b in zip(result, expected):
			assert np.array_equal(a, b)


@pytest.mark.parametrize("method", ["ipix2ang", "ipix2ang_center", "ipix2xy"])
def test_uint64_at_the_top_of_the_range(method):
	'''
	Test uint64 input at depth 30, where nbins is above the int64 range of an
	ordinary catalog column and the last valid ipix is 6*2^60-1.

	uint64 is the one integer dtype whose range does not fit in int64, so it is
	narrowed explicitly rather than left to NumPy's promotion rules (which changed
	in NumPy 2 and, if they landed on float64, would have a spacing of 1024 here).
	'''
	q = QLSC(depth=30)
	call = getattr(q, method)

	last = q.nbins - 1
	result = call(np.array([last], dtype=np.uint64))
	for a, b in zip(result, call(np.array([last], dtype=np.int64))):
		assert np.array_equal(a, b)

	# just past the end of the scheme, past the int64 range, and the largest uint64
	for bad in (q.nbins, q.nbins + 1, 2**63, 2**63 + 5, 2**64 - 1):
		with pytest.raises(ValueError):
			call(np.array([bad], dtype=np.uint64))


@pytest.mark.parametrize("method", ["ipix2ang", "ipix2ang_center", "ipix2xy"])
def test_non_integer_sequences_raise(method):
	'''
	Test that a Python sequence of non-integers is refused rather than converted.

	NumPy's "safe cast" rule refuses an existing float64 array but does not apply
	to an array it is constructing for the first time, so a list literal is the
	form in which [1.5, 2.5] could reach the pixel lookup as [1, 2].
	'''
	q = QLSC(depth=5)
	call = getattr(q, method)

	for bad in ([1.5, 2.5], (1.9, 2.9), [1.0, 2.0], [True, False], ["10", "20"],
				[None, None], [2**70], [[1, 2], [3]]):
		with pytest.raises((ValueError, TypeError)):
			call(bad)


@pytest.mark.parametrize("method", ["ipix2ang", "ipix2ang_center", "ipix2xy"])
def test_empty_input(method):
	'''
	Test that empty input returns empty arrays rather than raising.
	'''
	q = QLSC(depth=5)
	for form in ([], np.array([], dtype=np.int64), np.array([])):
		result = getattr(q, method)(form)
		for item in result:
			assert isinstance(item, np.ndarray)
			assert item.shape == (0,)


@pytest.mark.parametrize("method", ["ipix2ang", "ipix2ang_center", "ipix2xy"])
def test_out_of_range_values_raise(method):
	'''
	Test that out-of-range ipix values raise instead of decoding to coordinates
	that belong to no pixel of the scheme.
	'''
	q = QLSC(depth=5)
	call = getattr(q, method)

	for bad in (-1, q.nbins, q.nbins + 100):
		with pytest.raises(ValueError):
			call(bad)
		with pytest.raises(ValueError):
			call(np.array([0, bad]))

	# the message must name the offending value so it can be found in a big array
	with pytest.raises(ValueError, match=str(q.nbins)):
		call(np.array([0, 1, q.nbins, 2]))


@pytest.mark.parametrize("method", ["ipix2ang", "ipix2ang_center", "ipix2xy"])
def test_non_integer_values_raise(method):
	'''
	Test that float (and other non-integer) input raises rather than being
	silently truncated to a different pixel.
	'''
	q = QLSC(depth=5)
	call = getattr(q, method)

	for bad in (1.5, 10.0, np.float64(10.0), "10", np.array([1.5, 2.5]), np.array([1.0, 2.0]), True):
		with pytest.raises(ValueError):
			call(bad)


def test_c_level_rejects_float_array():
	'''
	Test that the C wrapper itself refuses a float array rather than truncating
	it (the Python layer rejects it first, so this guards the C contract).
	'''
	q = QLSC(depth=5)
	with pytest.raises(TypeError):
		q3c.ipix2ang(q._hprm, np.array([1.5, 2.5]))
	with pytest.raises(TypeError):
		q3c.ipix2xy(q._hprm, np.array([1.5, 2.5]))


@pytest.mark.parametrize("bad", [
	pytest.param([1.5, 2.5], id="list-of-floats"),
	pytest.param((1.9, 2.9), id="tuple-of-floats"),
	pytest.param([1.0, 2.0], id="list-of-whole-floats"),
	pytest.param([True, False], id="list-of-bools"),
	pytest.param(["10", "20"], id="list-of-strings"),
	pytest.param([None], id="list-of-none"),
	pytest.param([2**70], id="list-beyond-int64")
])
def test_c_level_rejects_non_integer_sequences(bad):
	'''
	Test that the C wrapper refuses a Python sequence of non-integers.

	This is the gap a plain "give me an int64 array" conversion leaves: NumPy's
	safe-cast rule refuses an existing float64 array, but it will happily BUILD an
	int64 array from a list of floats, booleans or numeric strings, so
	q3c.ipix2ang(hprm, [1.5, 2.5]) returned the centers of pixels 1 and 2 without
	a word. The dtype of the input is now settled before anything is converted.
	'''
	q = QLSC(depth=5)
	with pytest.raises(TypeError):
		q3c.ipix2ang(q._hprm, bad)
	with pytest.raises(TypeError):
		q3c.ipix2xy(q._hprm, bad)


def test_c_level_accepts_unsigned_arrays():
	'''
	Test that the C wrapper accepts unsigned integer arrays, uint64 included: an
	int64 conversion refuses uint64 as an unsafe cast, but every valid ipix fits.
	uint64 values that do not fit in int64 land in the negative half and are
	caught by the bounds check.
	'''
	q = QLSC(depth=5)
	expected = q3c.ipix2ang(q._hprm, np.array([3, 17, 4095], dtype=np.int64))

	for dtype in (np.uint16, np.uint32, np.uint64): # 4095 does not fit in uint8
		result = q3c.ipix2ang(q._hprm, np.array([3, 17, 4095], dtype=dtype))
		for a, b in zip(result, expected):
			assert np.array_equal(a, b)

	with pytest.raises(ValueError):
		q3c.ipix2ang(q._hprm, np.array([2**63 + 5], dtype=np.uint64))


def test_c_level_bounds_check():
	'''
	Test that the C wrapper bounds checks ipix values in both the scalar and
	array forms, independently of the Python-level check.
	'''
	q = QLSC(depth=5)
	for bad in (-1, q.nbins):
		with pytest.raises(ValueError):
			q3c.ipix2ang(q._hprm, bad)
		with pytest.raises(ValueError):
			q3c.ipix2xy(q._hprm, bad)
		with pytest.raises(ValueError):
			q3c.ipix2ang(q._hprm, np.array([0, bad]))
		with pytest.raises(ValueError):
			q3c.ipix2xy(q._hprm, np.array([0, bad]))


def test_large_array():
	'''
	Test a large array end to end: the whole point of vectorizing is bulk use.
	'''
	q = QLSC(depth=30)
	rng = np.random.default_rng(seed=20260830)
	ipix = rng.integers(0, q.nbins, size=100_000, dtype=np.int64)

	ra, dec = q.ipix2ang(ipix)
	assert ra.shape == dec.shape == ipix.shape
	assert np.isfinite(ra).all() and np.isfinite(dec).all()
	assert ((0 <= ra) & (ra <= 360)).all()
	assert ((-90 <= dec) & (dec <= 90)).all()
	assert np.array_equal(q.ang2ipix(ra=ra, dec=dec), ipix)

	facenum, x, y = q.ipix2xy(ipix)
	assert ((0 <= facenum) & (facenum <= 5)).all()
	assert ((-1 <= x) & (x <= 1)).all()
	assert ((-1 <= y) & (y <= 1)).all()


def test_ipix2polygon_remains_scalar():
	'''
	Test that ipix2polygon is unchanged (documented as scalar-only) and still
	rejects array input loudly.
	'''
	q = QLSC(depth=5)
	polygon = q.ipix2polygon(10)
	assert polygon.shape == (4, 2)

	with pytest.raises((ValueError, TypeError)):
		q.ipix2polygon(np.array([10, 11]))
