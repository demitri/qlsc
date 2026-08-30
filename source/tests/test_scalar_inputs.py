
import warnings

import pytest
import numpy as np

from qlsc import QLSC, QLSCIndex, distance, sindist, xy2facenum

# Every public API that takes scalar angle/coordinate values must accept
# scalar-LIKE values (NumPy scalars, 0-d arrays, size-1 arrays) without relying
# on NumPy's deprecated size-1-array-to-scalar conversion, which will become an
# error in a future NumPy release. This enumerates the class; a new scalar-taking
# API should be added to the list. (add_points is deliberately absent: it is an
# array API by contract and coerces its inputs with np.asarray itself.)

q30 = QLSC(depth=30)
POLYGON = np.array([[5.,-5.],[5.,5.],[-5.,5.]])

def _query_index(wrap):
	idx = QLSCIndex(qlsc=q30)
	idx.add_point(10., 30.)
	return idx.radial_query(ra=wrap(10.), dec=wrap(30.), radius=wrap(1.))

def _add_point_index(wrap):
	idx = QLSCIndex(qlsc=q30)
	idx.add_point(wrap(10.), wrap(30.))
	# Query the point back and return the matches: a mis-stored point (e.g. a NumPy
	# array bound by sqlite3 as a BLOB) still increments number_of_points, so the
	# row count is exactly the property that stays correct while the data is wrong.
	matches = idx.radial_query(ra=10., dec=30., radius=0.5)
	assert matches.shape == (2,), f"point added via add_point was not found by radial_query (matches shape: {matches.shape})"
	return matches

# name, callable taking a wrap function applied to every scalar argument
scalar_api_calls = [
	("ang2ipix", lambda w: q30.ang2ipix(w(10.), w(30.))),
	("ang2ipix_xy", lambda w: q30.ang2ipix_xy(w(10.), w(30.))),
	("ang2ipix_at_depth", lambda w: q30.ang2ipix_at_depth(w(10.), w(30.), 10)),
	("face_number", lambda w: q30.face_number(w(10.), w(30.))),
	("xy2ang", lambda w: q30.xy2ang(facenum=1, x=w(0.1), y=w(0.2))),
	("xy2ipix", lambda w: q30.xy2ipix(facenum=1, x=w(0.1), y=w(0.2))),
	("point_in_polygon", lambda w: q30.point_in_polygon(ra=w(0.), dec=w(0.), polygon=POLYGON)),
	("distance", lambda w: distance(w(0.), w(0.), w(30.), w(0.))),
	("sindist", lambda w: sindist(w(0.), w(0.), w(30.), w(0.))),
	# y beyond +1 crosses onto the top face, so the result (0) differs from the
	# input facenum - an implementation echoing its argument would be caught
	("xy2facenum", lambda w: xy2facenum(w(0.1), w(1.5), 1)),
	("radial_query", _query_index),
	("add_point", _add_point_index)
]

scalar_like_wraps = [
	pytest.param(lambda v: np.array([v]), id="size-1-array"),
	pytest.param(lambda v: np.array(v), id="0d-array"),
	pytest.param(np.float64, id="numpy-scalar")
]

@pytest.mark.parametrize("wrap", scalar_like_wraps)
@pytest.mark.parametrize("name, call", scalar_api_calls, ids=[n for n, c in scalar_api_calls])
def test_scalar_apis_accept_scalar_like(name, call, wrap):
	'''
	Test that scalar-taking APIs accept scalar-like values, give the same answer
	as with plain floats, and emit no deprecation warnings.
	'''
	with warnings.catch_warnings():
		warnings.simplefilter("error", DeprecationWarning)
		result = call(wrap)
	expected = call(lambda v: v)

	if isinstance(expected, dict):
		assert result == expected
	else:
		# some APIs legitimately return an array for array input and a scalar
		# for scalar input (e.g. ang2ipix), so compare values elementwise
		assert np.all(np.asarray(result) == np.asarray(expected))

def test_scalar_apis_reject_multi_element_arrays():
	'''
	Test that a scalar-only API rejects a multi-element array with a clear error.
	'''
	with pytest.raises(TypeError):
		q30.ang2ipix_xy(np.array([10., 20.]), np.array([30., 40.]))
	with pytest.raises(TypeError):
		distance(np.array([0., 1.]), 0., 30., 0.)
