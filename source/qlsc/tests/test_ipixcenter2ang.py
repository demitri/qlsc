
import pytest
from numpy.testing import assert_approx_equal

from qlsc import QLSC

# bin_level, ipix, ra, dec
expected_results_30 = [
	(2, 0, 315.0, 43.31385665828306),
	(2, 3, 315.0, 70.52877936550931),
	(2, 5, 45.0, 43.31385665828306),
	(2, 9, 225.0, 70.52877936550931),
	(2, 15, 135.0, 43.31385665828306),
	(2, 17, 345.96375653207355, -36.039893430303856),
	(2, 22, 14.036243467926479, -13.63302222536641),
	(2, 25, 345.96375653207355, 13.63302222536641),
	(2, 30, 14.036243467926479, 36.039893430303856),
	(2, 82, 251.56505117707798, -51.67118189854412),
	(2, 83, 225.0, -70.52877936550931),
	(2, 87, 108.43494882292202, -51.67118189854412),
	(2, 90, 315.0, -43.31385665828306),
	(2, 95, 45.0, -43.31385665828306)
]

expected_results = expected_results_30

@pytest.mark.parametrize("bin_level, ipix, ra, dec", expected_results)
def test_ipix2ang_center(bin_level, ipix, ra, dec):
	'''
	Test QLSC ipixcenter2ang.
	'''
	q = QLSC(bin_level=bin_level) # use default nside value
	ra_out, dec_out = q.ipix2ang_center(ipix)
	
	assert_approx_equal(ra, ra_out)
	assert_approx_equal(dec, dec_out)
