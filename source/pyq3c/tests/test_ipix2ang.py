
import pytest
from numpy.testing import assert_approx_equal

from pyq3c import Q3C

# ipix, ra, dec
expected_results = [
	(294541929365396512, 12.3, 45.6),
	(6438254810785524898, 311.76, -84.33),
	(2017612633061982208, 0, 0),
	(864691128455135232, 180, 90)
]

@pytest.mark.parametrize("ipix, ra, dec", expected_results)
def test_ipix2ang(ipix, ra, dec):
	'''
	Test Q3C ipix2ang.
	'''
	q3c = Q3C() # use default nside value
	ra_out, dec_out = q3c.ipix2ang(ipix)
	
	assert_approx_equal(ra, ra_out)
	assert_approx_equal(dec, dec_out)
