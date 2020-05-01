
import pytest
import numpy as np
from numpy.testing import assert_approx_equal # use foe scalars
from numpy.testing import assert_allclose     # use for arrays

from qlsc.utilities import _normalize_ang

ra_int = np.array([-720, -359, -45,  0, 45,  90, 135, 270, 359, 400, 720], dtype=np.int32)
ra = np.array([-720, -359, -45,  0, 45,  90, 135, 270, 359, 400, 720], dtype=np.double)
dec= np.array([ -45,   37,  45,-95, 65, 135,  54,  89, 9.5,  33,   0], dtype=np.double)

ra_norm  = np.array([  0.,   1., 315., 180.,  45., 270., 135., 270., 359.,  40.,  0.])
dec_norm = np.array([-45.,  37.,  45., -85.,  65.,  45.,  54.,  89.,   9.5, 33.,  0. ])

# ra, dec, norm_ra, norm_dec
scalar_expected_results = [
	(-734.,  -45, 346., -45.),
	( 175., 179., 355.,   1.),
	( 837.,-451., 297., -89.),
	( 145., -34., 145., -34.),
	( 200 ,  23 , 200 ,  23 )
]

def test_normalize_ang_individual_arrays():
	'''
	Test _normalize_ang when given individual ra, dec arrays.
	'''
	ra_id = id(ra)
	dec_id = id(dec)
	
	_normalize_ang(ra, dec)
	
	assert_allclose(ra, ra_norm, err_msg="The RA array did not normalize correctly.")
	assert_allclose(dec, dec_norm, err_msg="The dec array did not normalize correctly.")
	
	# the arrays should not have been copied
	assert ra_id == id(ra), "the ra array was unexpectedly copied"
	assert dec_id == id(dec), "the dec array was unexpectedly copied"

	# ultimately the point...	
	assert (ra >= 0).all() and (ra <= 360).all(), "ra values out of [0,360]!!"
	assert (dec >= -90).all() and (dec <= 90).all(), "dec values out of [-90,90]!!"

def test_normalize_ang_int_arrays():
	'''
	Test handling of non-floating point arrays.
	'''
	# ra is np.int32 and we asked not to copy the values
	with pytest.raises(Exception):
		_normalize_ang(ra_int, dec, copy=False)

def test_normalize_ang_points_array():
	'''
	Test passing in a 'points' array.
	'''
	
	# create a 2D array; only use the first and last points
	points = np.squeeze(np.dstack((np.copy(ra).astype(np.double), dec)))
	_normalize_ang(points=points)
	
	# change shape of output
	points_norm = np.squeeze(np.dstack((ra_norm, dec_norm)))
	
	assert_allclose(points[:,0], points_norm[:,0], err_msg="converting a 'points' array failed; ra didn't match")
	assert_allclose(points[:,1], points_norm[:,1], err_msg="converting a 'points' array failed; dec didn't match")

	assert (points[:,0] >= 0).all() and (points[:,0] <= 360).all(), "ra values out of [0,360]!!"
	assert (points[:,1] >= -90).all() and (points[:,1] <= 90).all(), "dec values out of [-90,90]!!"
	
def test_normalize_ang_parameters():
	'''
	Test that the correct parameters are used.
	'''
	
	with pytest.raises(Exception):
		_normalize_ang(ra=ra, points=points)

def test_normalize_ang_copy_made():
	'''
	Make sure that arrays are not modified if copies were requested.
	'''
	
	ra_id = id(ra)
	dec_id = id(dec)

@pytest.mark.parametrize("ra, dec, ra_norm, dec_norm", scalar_expected_results)
def test_normalize_ang_scalars(ra, dec, ra_norm, dec_norm):
	'''
	Test that scalar values work.
	'''
	ra_out, dec_out = _normalize_ang(ra, dec)
	assert_approx_equal(ra_out, ra_norm, err_msg="RA scalar values did not match.")
	assert_approx_equal(dec_out, dec_norm, err_msg="dec scalar values did not match.")

def test_normalize_ang_float_arrays():
	'''
	Test float arrays when copied.
	'''
	ra = np.array([451,23,54,-5], dtype=np.float32)
	dec = np.array([12,54,-95,5], dtype=np.float32)
	
	ra_in_id =  id(ra)
	dec_in_id = id(dec)
	
	ra_out, dec_out = _normalize_ang(ra, dec, copy=True)
	
	assert_allclose(ra_out, np.array([91., 23., 234., 355.], dtype=np.float32), err_msg="ra doesn't match when input is float array and copy requested")
	assert_allclose(dec_out, np.array([12., 54., -85., 5.], dtype=np.float32), err_msg="dec doesn't match when input is float array and copy requested")
	
	assert ra_in_id != id(ra_out), "ra copy requested, but original modified"
	assert dec_in_id != id(dec_out), "dec copy requested, but original modified"
	
def test_normalize_ang_points_return_on_copy():
	'''
	Test that a single 'points' array is returned when one points array is provided.
	'''
	points = np.squeeze(np.dstack((np.array([451,23,54,-5], dtype=np.float32), np.array([12,54,-95,5], dtype=np.float32))))
	
	points_id = id(points)
	
	points_out = _normalize_ang(points=points, copy=True)
	
	points_norm = np.array([[ 91.,  12.],
							[ 23.,  54.],
							[234., -85.],
							[355.,   5.]], dtype=np.float32)
	
	assert_allclose(points_out, points_norm, err_msg="points copy test doesn't match")
	
	assert points_id != id(points_out), "the original 'points' array was modified, not copied"
	
def test_normalize_ang_mismatching_array_lengths():
	'''
	Test mismatching array lengths.
	'''
	ra = np.array([1,2,3,4], dtype=np.double)
	dec = np.array([53,654], dtype=np.double)
	
	with pytest.raises(Exception):
		_normalize_ang(ra,dec)
	
	
	
	