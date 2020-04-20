
import pytest
from numpy.testing import assert_approx_equal

from qlsc import QLSC

# bin_level, ipix, ra, dec
expected_results_30 = [
	(30, 294541929365396512, 12.3, 45.6),
	(30, 6438254810785524898, 311.76, -84.33),
	(30, 2017612633061982208, 0, 0),
	(30, 864691128455135232, 180, 90)
]

expected_results = expected_results_30

@pytest.mark.parametrize("bin_level, ipix, ra, dec", expected_results)
def test_ipix2ang(bin_level, ipix, ra, dec):
	'''
	Test QLSC ipix2ang.
	'''
	q = QLSC(bin_level=bin_level) # use default nside value
	ra_out, dec_out = q.ipix2ang(ipix)
	
	assert_approx_equal(ra, ra_out)
	assert_approx_equal(dec, dec_out)
