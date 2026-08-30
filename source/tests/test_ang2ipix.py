
import pytest
import numpy as np
#from numpy.testing import assert_approx_equal # use foe scalars
#from numpy.testing import assert_allclose     # use for arrays

from qlsc import QLSC

# depth, ra, dec, ipix
expected_results = [
	(0, 0, 0, 1),
	(0, -22, 0, 1),
	(0, -45, 0, 1),
	(0, 27, 0, 1),
	(0, 45, 0, 2),
	(0, 46, 0, 2),
	(0, 256, 0, 4),
	(0, 12, 34, 1),
	(0, 12, 87, 0),
	(0, 12, -73, 5),
	(0, 12, -3, 1),
	(30, 0, 0, 2017612633061982208),
	(30, 180, 0, 4323455642275676160),
	(30, -1, 0, 1825388318008017157),
	(30, 359, 0, 1825388318008017157),
	(30, 0, -89.9, 6629299373591235202),
	(30, 57.3, 88, 480333897118235422)
]

# depth, ra, dec, ipix
expected_results_q0 = [
	(0, np.array([0, -22, -45, 27, 45, 46, 256, 12, 12, 12, 12], dtype=np.double), np.array([0, 0, 0, 0, 0, 0, 0, 34, 87, -73, -3], dtype=np.double), np.array([1, 1, 1, 1, 2, 2, 4, 1, 0, 5, 1], dtype=np.int64)),
	(0, np.array([0, -22, -45, 27, 45, 46, 256, 12, 12, 12, 12], dtype=np.double), np.array([0, 0, 0, 0, 0, 0, 0, 34, 87, -73, -3], dtype=np.double), np.array([1, 1, 1, 1, 2, 2, 4, 1, 0, 5, 1])),
	(30, np.array([0., 180., -1., 359., 0., 57.3], dtype=np.double), np.array([0., 0., 0., 0., -90., 88.], dtype=np.double), np.array([2017612633061982208, 4323455642275676160, 1825388318008017157, 1825388318008017157, 6629298651489370112, 480333897118235422], dtype=np.int64)),
	(30, np.array([0., 180., -1., 359., 0., 57.3], dtype=np.double), np.array([0., 0., 0., 0., -90., 88.], dtype=np.double), np.array([2017612633061982208, 4323455642275676160, 1825388318008017157, 1825388318008017157, 6629298651489370112, 480333897118235422], dtype=np.int64))
]

@pytest.mark.parametrize("depth, ra, dec, ipix", expected_results)
def test_ang2ipix_scalar(depth, ra, dec, ipix):
	'''
	Test QLSC ang2ipix using scalar values.
	'''
	q = QLSC(depth=depth)
	assert ipix == q.ang2ipix(ra, dec)

@pytest.mark.parametrize("depth, ra, dec, ipix", expected_results_q0)
def test_ang2ipix_array(depth, ra, dec, ipix):
	'''
	Test ang2ipix using arrays.
	'''
	q = QLSC(depth=depth)
	assert np.array_equal(ipix, q.ang2ipix(ra, dec))

# Only the depth=0 rows are used here: their ra,dec values are integers, which are exactly
# representable in float32; a value like 57.3 changes when cast, legitimately changing the
# expected ipix at high depths.
@pytest.mark.parametrize("depth, ra, dec, ipix", [r for r in expected_results_q0 if r[0] == 0])
def test_ang2ipix_array_float32(depth, ra, dec, ipix):
	'''
	Test ang2ipix using float32 arrays, which the method must convert to double.
	'''
	q = QLSC(depth=depth)
	assert np.array_equal(ipix, q.ang2ipix(ra.astype(np.float32), dec.astype(np.float32)))

def test_ang2ipix_0d_array():
	'''
	Test that 0-d array inputs work (regression: the length-1-array guard
	indexed .shape[0] unconditionally, raising IndexError for 0-d arrays).
	'''
	q = QLSC(depth=30)
	assert q.ang2ipix(np.array(10.), np.array(30.)) == q.ang2ipix(10., 30.)

def test_ang2ipix_dec_below_range():
	'''
	Test that out of range dec values (< -90) are caught.
	'''
	# The underlying Q3C code will truncate dec values outside of [-90,90].
	# Make sure that we check for this before calling the function.
	ra = 56.0
	dec = -135
	
	q = QLSC(depth=30)
	
	truncated_dec_ipix = q.ang2ipix(56.0, -90) # 6629298651489370112
	correct_ipix = 5816706551187910542
	
	assert q.ang2ipix(ra,dec) == correct_ipix, "QLSC.ang2ipix is allowing out of range dec values (< -90) to be truncated."
	
def test_ang2ipix_dec_above_range():
	'''
	Test that out of range dec values (> 90) are caught.
	'''
	# The underlying Q3C code will truncate dec values outside of [-90,90].
	# Make sure that we check for this before calling the function.
	ra = 56.0
	dec = 255.0
	
	q = QLSC(depth=30)
	
	truncated_dec_ipix = q.ang2ipix(56.0, 90) # 864691128455135232
	correct_ipix = 6037741812941379589
	
	assert q.ang2ipix(ra,dec) == correct_ipix, "QLSC.ang2ipix is allowing out of range dec values (> 90) to be truncated."

def test_ang2ipix_points():
	'''
	Test ang2ipix from array.
	'''
	points = np.array([[0, 12], [35, 32], [32, -32]], dtype=np.double)

	q = QLSC(depth=30)
	
	ipix = q.ang2ipix(points=points)
	
	assert (ipix == np.array([2029048250313091210, 2275510515682463237, 1550784043278934694], dtype=np.int64)).all(), "use of points array failed"
	
def test_ang2ipix_points_with_dec_out_of_range():
	'''
	Test ang2ipix from array with dec values out of range.
	'''
	points = np.array([[0, 102], [35, -132], [32, 99]], dtype=np.double)

	q = QLSC(depth=30)
	
	ipix = q.ang2ipix(points=points)
	
	assert (ipix == np.array([ 876126745706244234, 5824629191300107707,  680379582899038849], dtype=np.int64)).all(), "use of points array with dec out of range failed"

