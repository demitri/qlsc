
import pytest

from pyq3c import Q3C

# nside, ra, dec, ipix
expected_results = [
	(1073741824, 0, 0, 2017612633061982208),
	(1073741824, 0, 180, 480383960252852906),
	(1073741824, -1, 0, 1825388318008017157),
	(1073741824, 359, 0, 1825388318008017157),
	(1073741824, 0, -90, 6629298651489370112),
]

@pytest.mark.parametrize("nside, ra, dec, ipix", expected_results)
def test_ang2ipix(nside, ra, dec, ipix):
	'''
	Test Q3C ang2ipix.
	'''
	q3c = Q3C(nside=nside)
	assert ipix == q3c.ang2ipix(ra, dec)
